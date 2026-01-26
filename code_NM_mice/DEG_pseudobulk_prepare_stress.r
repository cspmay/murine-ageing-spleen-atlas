library(tidyverse)
library(scran)
library(Matrix)
library(glmGamPoi)

out_dir <- './NM_mice/out'

DE_out_dir <- file.path(out_dir, 'pseudobulk_DE')
dir.create(DE_out_dir, showWarnings = FALSE)

cell_metadata <- file.path(out_dir, 'annotated_metadata_NM_mice.tsv') |>
# cell_metadata <- file.path(out_dir, 'annotated_T_metadata_NM_mice.tsv') |>
# cell_metadata <- file.path(out_dir, 'annotated_B_metadata_NM_mice.tsv') |>
# cell_metadata <- file.path(out_dir, 'annotated_Gzmk_metadata_NM_mice.tsv') |>
  read.csv(sep = '\t', row.names = 1) |>
  select(sample_name, age_group, condition, cell_type_V2, dataset) |>
  mutate(condition = factor(condition, levels = c('control', 'stress')))
cell_metadata |> head(3)

all_cell_types <- unique(cell_metadata$cell_type_V2)
all_cell_types

all_conditions <- unique(cell_metadata$condition)
all_conditions

counts <- file.path(out_dir, 'all_annotated', 'matrix.mtx') |>
# counts <- file.path(out_dir, 'T_annotated', 'matrix.mtx') |>
# counts <- file.path(out_dir, 'B_annotated', 'matrix.mtx') |>
# counts <- file.path(out_dir, 'Gzmk_annotated', 'matrix.mtx') |>
  readMM() |>
  t()
features <- file.path(out_dir, 'all_annotated', 'features.txt') |>
# features <- file.path(out_dir, 'T_annotated', 'features.txt') |>
# features <- file.path(out_dir, 'B_annotated', 'features.txt') |>
# features <- file.path(out_dir, 'Gzmk_annotated', 'features.txt') |>
  readLines()
rownames(counts) <- features
barcodes <- file.path(out_dir, 'all_annotated', 'barcodes.txt') |>
# barcodes <- file.path(out_dir, 'T_annotated', 'barcodes.txt') |>
# barcodes <- file.path(out_dir, 'B_annotated', 'barcodes.txt') |>
# barcodes <- file.path(out_dir, 'Gzmk_annotated', 'barcodes.txt') |>
  readLines()
colnames(counts) <- barcodes

counts <- as(counts, 'CsparseMatrix')

dim(counts)
dim(cell_metadata)
sum(rownames(cell_metadata) %in% colnames(counts))

cell_metadata <- cell_metadata |> 
  filter(age_group %in% c("old"))
sce <- SingleCellExperiment(
  assays = list(counts = counts[, rownames(cell_metadata)]),
  colData = cell_metadata
)
sce

export_pseudobulk <- function(sce, DE_out_dir) {
  print(paste0('Total cells: ', ncol(sce)))     
  cell_types_to_exclude <- c()
  excluded_cell_type_samples <- list()
  for(cell_type in all_cell_types) {  
    sample_cell_type_counts <- sce |>
      colData() |>
      as.data.frame() |>
      dplyr::count(sample_name, condition, cell_type_V2) |>
      filter(cell_type_V2 == cell_type)
    sufficient_groups <- sample_cell_type_counts |>
      filter(n >= 10) |>
      dplyr::count(condition) |>
      filter(n >= 3) |> 
      nrow()
    if (sufficient_groups != 2) {
      cell_types_to_exclude <- c(cell_types_to_exclude, cell_type)
    }
    samples_to_exclude <- sample_cell_type_counts |>
      filter(n < 10) |>
      pull(sample_name)        
    if (length(samples_to_exclude) > 0) {
      excluded_cell_type_samples[[cell_type]] <- samples_to_exclude 
    }
  }
  print(paste0('Excluded cell types: ', paste(cell_types_to_exclude, collapse = ', ')))
  filtered_sce <- sce[, !(colData(sce)$cell_type_V2 %in% cell_types_to_exclude)]  
  pseudobulk_sce <- filtered_sce |>
    pseudobulk(group_by = vars(sample_name, condition = condition, cell_type_V2))
  for(cell_type in names(excluded_cell_type_samples)) {
    samples_to_exclude <- excluded_cell_type_samples[[cell_type]]
    print(cell_type)
    print(paste0('+++ Excluded samples: ', paste(samples_to_exclude, collapse = ', ')))    
    pseudobulk_profiles_to_exclude <- !(
      (colData(pseudobulk_sce)$cell_type_V2 == cell_type) & 
      (colData(pseudobulk_sce)$sample_name %in% samples_to_exclude)
    )
    pseudobulk_sce <- pseudobulk_sce[, pseudobulk_profiles_to_exclude]
  }
  saveRDS(pseudobulk_sce, file.path(DE_out_dir, paste0('NM_mice_all_stress_pseudobulk.rds')))
  # saveRDS(pseudobulk_sce, file.path(DE_out_dir, paste0('NM_mice_T_stress_pseudobulk.rds')))
  # saveRDS(pseudobulk_sce, file.path(DE_out_dir, paste0('NM_mice_B_stress_pseudobulk.rds')))
  # saveRDS(pseudobulk_sce, file.path(DE_out_dir, paste0('NM_mice_Gzmk_stress_pseudobulk.rds')))
}
export_pseudobulk(sce, DE_out_dir)
