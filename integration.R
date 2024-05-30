suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(qs)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(future)))
suppressMessages(suppressWarnings(library(future.apply)))
suppressMessages(suppressWarnings(library(harmony)))
suppressMessages(suppressWarnings(library(Azimuth)))
options(future.globals.maxSize = 1000 * 64536^2)
plan(multisession, workers = 16, gc = TRUE)

process_seurat_object <- function(meta_path, input_path) {
  # Attempt to read and process the Seurat object
  tryCatch({
    cat("Loading and processing the Seurat object...\n")
    obj <- readRDS(input_path)

    # Join Layers and manipulate metadata
    obj <- JoinLayers(obj)
    obj_meta <- obj@meta.data %>%
      mutate(Sample = gsub("-Dejager", "_Dejager", orig.ident)) %>%
      rownames_to_column("Cell") %>%
      mutate(
        Barcode = str_extract(Cell, "[ATGC]{16}"),
        Clean_Sample = if_else(!is.na(Barcode), str_extract(Cell, paste0(".*", Barcode)), Sample),
        Sample_barcode = if_else(!is.na(Barcode), Clean_Sample, Cell)
      ) %>%
      group_by(Sample_barcode) %>%
      mutate(id = row_number(), Sample_barcode = if(n() > 1) paste(Sample_barcode, id, sep = "-") else Sample_barcode) %>%
      ungroup() %>%
      select(-id) %>%
      dplyr::select(Cell, Sample_barcode)

    # Metadata filtering and merging
    df_meta <- read.csv(meta_path)  # Assuming meta_path is a path to a CSV file
    filtered_metadata <- merge(obj_meta, df_meta %>% rename(Sample_barcode = Cell), by = "Sample_barcode") %>%
      column_to_rownames("Cell") %>%
      filter(!is.na(Batch) & !is.na(mitoRatio) & Batch %in% c(1, 2) & !is.na(Cohort)) %>%
      rename(orig_predicted.cluster = predicted.cluster, orig_predicted.cluster.score = predicted.cluster.score)

    obj <- subset(obj, cells = rownames(filtered_metadata))
    filtered_metadata <- filtered_metadata[rownames(obj@meta.data), ]
    obj@meta.data <- filtered_metadata

    # Seurat processing with Harmony
    obj <- obj %>%
      NormalizeData() %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
      ScaleData() %>%
      SCTransform(vars.to.regress = c("mitoRatio")) %>%
      RunPCA(assay = "SCT", npcs = 50)

    obj <- RunHarmony(obj, ncores = 16, 
                      group.by.vars = c("Sample", "Cohort", "Batch"), 
                      reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
    obj <- RunUMAP(obj, reduction = "harmony", assay = "SCT", dims = 1:50)
    obj <- FindNeighbors(object = obj, reduction = "harmony")
    obj <- FindClusters(obj, resolution = c(0.2, 0.4, 0.6, 0.8))

    # Running Azimuth and finding markers
    obj <- RunAzimuth(obj, reference = "humancortexref", assay = "SCT")
    sample.markers <- FindAllMarkers(obj, assay = "SCT", logfc.threshold = 0.1, test.use = "MAST", slot = "SCT", random.seed = 1, latent.vars = c("Sample"))

    # Save QS files
    qs::qsave(sample.markers, "sample.markers.qs")
    qs::qsave(obj, "harmonized_annot.qs")

    # Output paths of saved files
    list(Sample_Markers_Path = "sample.markers.qs", Harmonized_Annot_Path = "harmonized_annot.qs")
    
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
    NULL
  })
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("No file path provided. Please specify the path to a Seurat object file.")
} else {
  # Call the function with the provided file path
   meta_path <- args[1]  # First argument: Path to metadata CSV file
  input_path <- args[2]  # Second argument: Path to the input Seurat object file
   process_seurat_object(meta_path, input_path)
  
 
}
