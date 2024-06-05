process_seurat_object <- function(input_path, ref_path) {
  # Attempt to read and process the Seurat object
  tryCatch({
    cat("Loading and processing the Seurat object...\n")
    obj <- JoinLayers(qs::qread(input_path)) 
    cells_filtered <- rownames(obj@meta.data %>% filter(
      between(nCount_RNA, quantile(nCount_RNA, 0.25) - 3 * IQR(nCount_RNA), quantile(nCount_RNA, 0.75) + 3 * IQR(nCount_RNA)),
      between(nFeature_RNA, quantile(nFeature_RNA, 0.25) - 3 * IQR(nFeature_RNA), quantile(nFeature_RNA, 0.75) + 3 * IQR(nFeature_RNA)),
      between(mitoRatio, quantile(mitoRatio, 0.25) - 3 * IQR(mitoRatio), quantile(mitoRatio, 0.75) + 3 * IQR(mitoRatio))
    )))
    obj <- subset(obj, cells = cells_filtered)

    # Pre-process Seurat object with standard Seurat workflow
    obj.sample <- NormalizeData(obj)
    obj.sample <- FindVariableFeatures(obj.sample)
    obj.sample <- ScaleData(obj.sample)
    obj.sample <- RunPCA(obj.sample, nfeatures.print = 10)

    # Find significant PCs
    stdv <- obj.sample[["pca"]]@stdev
    sum.stdv <- sum(obj.sample[["pca"]]@stdev)
    percent.stdv <- (stdv / sum.stdv) * 100
    cumulative <- cumsum(percent.stdv)
    co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
    co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = T)[1] + 1
    min.pc <- min(co1, co2)

    # Finish pre-processing
    obj.sample <- RunUMAP(obj.sample, dims = 1:min.pc)
    obj.sample <- FindNeighbors(object = obj.sample, dims = 1:min.pc)
    obj.sample <- FindClusters(object = obj.sample, resolution = 0.2)
    annotations <- obj.sample@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp.poi <- (0.066 + 0.000757 * nrow(obj.sample@meta.data)) / 100 * nrow(obj.sample@meta.data)
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
    doublet.obj <- doubletFinder(seu = obj.sample, PCs = 1:min.pc, pK = (0.066 + 0.000757 * nrow(obj.sample@meta.data)) / 100, nExp = nExp.poi.adj, sct = TRUE)
    colnames(doublet.obj@meta.data)[8] <- "doublet_finder"
    obj.singlets <- subset(doublet.obj, doublet_finder == "Singlet")

    newobj <- obj.singlets %>%
      SCTransform(variable.features.n = 2000, conserve.memory = TRUE, vars.to.regress = c("mitoRatio")) %>%
      RunPCA(dims = 1:min.pc) %>%
      FindNeighbors(dims = 1:min.pc) %>%
      FindClusters(resolution = 0.2, algorithm = 1) %>%
      RunUMAP(dims = 1:min.pc)
    
    obj <- RunAzimuth(newobj, reference = "humancortexref", assay = "SCT")
    reference <- readRDS(ref_path)
    outputname <- unique(obj@meta.data$orig.ident)

    p1 <- DimPlot(reference, reduction = "refUMAP", group.by = "subclass", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
    p2 <- DimPlot(obj, reduction = "umap", group.by = "predicted.subclass", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle(paste(outputname, "Query transferred labels"))
    p <- p1 + p2

    metadata_csv <- obj@meta.data %>% rownames_to_column("Cell")

    # Save UMAP plot as PDF
    ggsave(paste0(outputname, "_UMAP.pdf"), p, width = 10, height = 8)
    save_path <- paste0(outputname, "_metadata.csv")
    write.csv(metadata_csv, save_path, quote = FALSE, row.names = FALSE)

    # Return the plots and save path of the RDS file
    list(UMAP = p1, RDS_Path = save_path)

  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
    NULL
  })
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No file path provided. Please specify the path to a Seurat object file.")
} else {
  # Call the function with the provided file path
  results <- process_seurat_object(args[1], args[2])
  print(results$UMAP)
}

