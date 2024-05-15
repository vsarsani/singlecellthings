suppressWarnings(library(Seurat))
suppressWarnings(library(patchwork))
suppressWarnings(library(sctransform))
suppressWarnings(library(qs))
#suppressWarnings(library(tidyverse))
#suppressWarnings(library(Azimuth))

process_seurat_object <- function(input_path) {
  # Attempt to read and process the Seurat object
  tryCatch({
    cat("Loading and processing the Seurat object...\n")
    obj <- qs::qread(input_path) %>%
      SCTransform(variable.features.n = 2000, conserve.memory = TRUE,
                  vars.to.regress = c("mitoRatio")) %>%
      RunPCA() %>%
      FindNeighbors(dims = 1:30) %>%
      FindClusters(resolution = 0.4, algorithm = 1) %>%
      RunUMAP(dims = 1:30)

    # Generate output name based on unique identifiers in metadata
    outputname <- unique(obj@meta.data$orig.ident)

    cat("Running Azimuth mapping...\n")
    obj <- RunAzimuth(obj, reference = "humancortexref", assay = "SCT")

    cat("Finding markers...\n")
    sample.markers <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE,
                                     test.use = "negbinom", recorrect_umi = FALSE)

    cat("Generating UMAP plot...\n")
    p1 <- DimPlot(obj, group.by = "predicted.subclass", label = TRUE, label.size = 3) +
      ggtitle(paste(outputname, "UMAP"))

    # Save UMAP plot as PDF
    ggsave(paste0(outputname, "_UMAP.pdf"), p1, width = 10, height = 8)

    cat("Filtering and preparing heatmap data...\n")
    sample.markers %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 1) %>%
      slice_head(n = 10) %>%
      ungroup() -> top10

    cat("Generating heatmap...\n")
    p2 <- DoHeatmap(obj, features = top10$gene, group.by = "predicted.subclass",
                    label = TRUE, size = 3) +
      ggtitle(paste(outputname, "Marker genes by class"))

    # Save heatmap plot as PDF
    ggsave(paste0(outputname, "_markerheatmap.pdf"), p2, width = 10, height = 8)

    # Save the processed Seurat object
    save_path = paste0(outputname, ".Rds")
    cat("Saving processed Seurat object to ", save_path, "\n")
    saveRDS(obj, file = save_path)

    # Return the plots and save path of the RDS file
    list(UMAP = p1, Heatmap = p2, RDS_Path = save_path)

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
  results <- process_seurat_object(args[1])
  print(results$UMAP)
  print(results$Heatmap)
}
