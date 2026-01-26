library(tidyverse)

out_dir = './out/scanpy_DE/'

output_dir <- "./out/common_genes_4ofall"
dir.create(output_dir, showWarnings = FALSE)

sig_files <- list.files(out_dir, pattern = "_sig.tsv$", recursive = TRUE, full.names = TRUE)

file_info <- tibble(
  file_path = sig_files,
  dataset = basename(dirname(sig_files)),
  cell_type = basename(sig_files) %>% sub("_sig.tsv", "", .)
)

cell_types <- unique(file_info$cell_type)

results_list <- map(cell_types, function(ct) {
  ct_files <- file_info %>% filter(cell_type == ct)
  genes_by_dataset <- map(ct_files$file_path, function(f) {
    df <- read_tsv(f, show_col_types = FALSE)
    up <- df %>% filter(lfc > 0) %>% pull(name)
    down <- df %>% filter(lfc < 0) %>% pull(name)
    list(up = up, down = down)
  })
  names(genes_by_dataset) <- ct_files$dataset
  all_up_genes <- map(genes_by_dataset, "up") %>% unlist()
  all_down_genes <- map(genes_by_dataset, "down") %>% unlist()
  up_counts <- table(all_up_genes)
  down_counts <- table(all_down_genes)
  min_datasets <- 4
  common_up <- names(up_counts[up_counts >= min_datasets])
  common_down <- names(down_counts[down_counts >= min_datasets])
  list(
    cell_type = ct,
    n_datasets = nrow(ct_files),
    common_up = common_up,
    common_down = common_down,
    n_common_up = length(common_up),
    n_common_down = length(common_down)
  )
})

names(results_list) <- cell_types

summary_df <- map_df(results_list, function(res) {
  tibble(
    cell_type = res$cell_type,
    n_datasets = res$n_datasets,
    n_common_up = res$n_common_up,
    n_common_down = res$n_common_down
  )
})
write_tsv(summary_df, file.path(output_dir, "summary_common_genes.tsv"))

all_common_up <- map_df(results_list, function(res) {
  if (length(res$common_up) > 0) {
    tibble(
      cell_type = res$cell_type,
      gene_name = res$common_up
    )
  } else {
    tibble(cell_type = character(), gene_name = character())
  }
})
all_common_down <- map_df(results_list, function(res) {
  if (length(res$common_down) > 0) {
    tibble(
      cell_type = res$cell_type,
      gene_name = res$common_down
    )
  } else {
    tibble(cell_type = character(), gene_name = character())
  }
})
write_tsv(all_common_up, 
          file.path(output_dir, "all_celltypes_common_upreg.tsv"))

write_tsv(all_common_down, 
          file.path(output_dir, "all_celltypes_common_downreg.tsv"))
