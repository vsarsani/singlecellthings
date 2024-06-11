suppressWarnings(library(SeuratDisk))
suppressWarnings(library(Seurat))
suppressWarnings(library(future))
suppressWarnings(library(qs))
suppressWarnings(library(tidyverse))
suppressWarnings(library(Azimuth))
suppressWarnings(library(DoubletFinder))

# Define the function for processing and classifying the Seurat object
process_and_classify <- function(input_path, ref_rds_path, ref_annoy_path, output_path) {
  tryCatch({
    cat("Loading and processing the Seurat object...\n")
    obj <- JoinLayers(qs::qread(input_path))
    
    # Filter cells based on IQR bounds
    cells_filtered <- rownames(obj@meta.data %>% filter(
      between(nCount_RNA, quantile(nCount_RNA, 0.25) - 3 * IQR(nCount_RNA), quantile(nCount_RNA, 0.75) + 3 * IQR(nCount_RNA)),
      between(nFeature_RNA, quantile(nFeature_RNA, 0.25) - 3 * IQR(nFeature_RNA), quantile(nFeature_RNA, 0.75) + 3 * IQR(nFeature_RNA)),
      between(mitoRatio, quantile(mitoRatio, 0.25) - 3 * IQR(mitoRatio), quantile(mitoRatio, 0.75) + 3 * IQR(mitoRatio))
    ))
    
    obj <- subset(obj, cells = cells_filtered)
    
    # Pre-process Seurat object with standard Seurat workflow
    obj.sample <- NormalizeData(obj) %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA(nfeatures.print = 10)
    
    # Find significant PCs
    stdv <- obj.sample[["pca"]]@stdev
    sum.stdv <- sum(obj.sample[["pca"]]@stdev)
    percent.stdv <- (stdv / sum.stdv) * 100
    cumulative <- cumsum(percent.stdv)
    co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
    co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = T)[1] + 1
    min.pc <- min(co1, co2)
    
    # Finish pre-processing
    obj.sample <- RunUMAP(obj.sample, dims = 1:min.pc) %>%
      FindNeighbors(dims = 1:min.pc) %>%
      FindClusters(resolution = 0.2)
    
    annotations <- obj.sample@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp.poi <- (0.066 + 0.000757 * nrow(obj.sample@meta.data)) / 100 * nrow(obj.sample@meta.data)
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
    doublet.obj <- doubletFinder(seu = obj.sample, PCs = 1:min.pc, pK = (0.066 + 0.000757 * nrow(obj.sample@meta.data)) / 100, nExp = nExp.poi.adj, sct = TRUE)
    colnames(doublet.obj@meta.data)[8] <- "doublet_finder"
    obj.singlets <- subset(doublet.obj, doublet_finder == "Singlet")
    
    newobj <- obj.singlets %>%
      SCTransform(variable.features.n = 2000, ncells=6000,conserve.memory = TRUE, vars.to.regress = c("mitoRatio")) %>%
      RunPCA(dims = 1:30) %>%
      FindNeighbors(dims = 1:30) %>%
      FindClusters(resolution = 0.2, algorithm = 1) %>%
      RunUMAP(dims = 1:30)
    
    obj <- newobj
    
    # Continue with classification
    cat("Loading the reference map...\n")
    ref.names <- list(
    map = ref_rds_path,
    ann = ref_annoy_path
  )
    map <- readRDS(ref.names$map)
    map[["refdr.annoy.neighbors"]] <- LoadAnnoyIndex(
    object = map[["refdr.annoy.neighbors"]],
    file = file.path(ref.names$ann)
  )

    reference <- map
    dims <- as.double(slot(reference, "neighbors")$refdr.annoy.neighbors@alg.info$ndim)
    meta.data <- names(slot(reference, "meta.data"))
    NNTransform <- function(
  object,
  meta.data,
  neighbor.slot = "query_ref.nn",
  key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}



    
    
    anchors <- FindTransferAnchors(
      reference = reference,
      query = obj,
      k.filter = NA,
      reference.neighbors =  "refdr.annoy.neighbors",
      reference.assay = "refAssay",
      query.assay = "SCT",
      reference.reduction = "refDR",
      normalization.method = "SCT",
      features = rownames(Loadings(reference[["refDR"]])),
      dims = 1:dims,
      n.trees = 20,
      mapping.score.k = 100,
      verbose = TRUE
    )
    
    annotation.levels <- c('Subclass', 'Class_8', 'Supertype')
    refdata <- lapply(X = annotation.levels, function(x) {
      reference[[x, drop = TRUE]]
    })
    names(refdata) <- annotation.levels
    
    query <- TransferData(
      reference = reference,
      query = obj,
      query.assay = "SCT",
      dims = 1:dims,
      anchorset = anchors,
      refdata = refdata,
      n.trees = 20,
      store.weights = TRUE,
      k.weight = 50,
      verbose = TRUE
    )
    
    query <- IntegrateEmbeddings(
      anchorset = anchors,
      reference = reference,
      query = query,
      query.assay = "SCT",
      reductions = "pcaproject",
      reuse.weights.matrix = TRUE,
      verbose = TRUE
    )
    
    query[["query_ref.nn"]] <- FindNeighbors(object = Embeddings(reference[["refDR"]])[, 1:dims], query = Embeddings(query[["integrated_dr"]]), return.neighbor = TRUE, l2.norm = TRUE, verbose = TRUE)
    
    query <- NNTransform(
      object = query,
      meta.data = reference[[]]
    )
    
    # Project the query to the reference UMAP
    query[["ref.umap"]] <- RunUMAP(
      object = query[["query_ref.nn"]],
      reduction.model = reference[["refUMAP"]],
      reduction.key = 'UMAP_',
      verbose = TRUE, assay = "SCT"
    )
    
    query <- AddMetaData(
      object = query,
      metadata = MappingScore(anchors = anchors, ndim = dims),
      col.name = "mapping.score"
    )
    
    # Save the classified Seurat object
    qsave(query, output_path)
    cat("Classified Seurat object saved to", output_path, "\n")
    
    return(query)
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
    NULL
  })
}

# Main script to parse arguments and call the function
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Please provide the input Seurat object path, the reference RDS path, the reference Annoy path, and the output path.")
} else {
  input_path <- args[1]
  ref_rds_path <- args[2]
  ref_annoy_path <- args[3]
  output_path <- args[4]
  result <- process_and_classify(input_path, ref_rds_path, ref_annoy_path, output_path)
}
