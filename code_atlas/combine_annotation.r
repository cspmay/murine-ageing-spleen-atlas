library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

out_dir <- './atlas/out/'

V2_metadata <- file.path(out_dir, 'annotated_metadata.tsv') |>
    read_delim(show_col_types = F)  |>
    rename(barcode = `...1`) |>
    select(barcode, sample_name, age_group, condition, dataset, cell_type_V2)
V2_metadata |> head(3)

unique(V2_metadata$cell_type_V2)

B_V2_metadata <- file.path(out_dir, 'annotated_B_metadata.tsv') |>
	read_delim(show_col_types = F) |>
    rename(barcode = `...1`) |>
    select(barcode, sample_name, age_group, condition, dataset, cell_type_V2)
B_V2_metadata |> head(3)

unique(B_V2_metadata$cell_type_V2)

T_V2_metadata <- file.path(out_dir, 'annotated_T_metadata.tsv') |>
	read_delim(show_col_types = F)  |>
    rename(barcode = `...1`) |>
    select(barcode, sample_name, age_group, condition, dataset, cell_type_V2)
T_V2_metadata |> head(3)

all_except_subsets_metadata <- V2_metadata |>
  filter(!(barcode %in% T_V2_metadata$barcode)) |>
  filter(!(barcode %in% B_V2_metadata$barcode))

nrow(all_except_subsets_metadata)

T_contaminating_barcodes <- file.path(out_dir, 'T_contaminating_barcodes.txt') |>
  readLines()
length(T_contaminating_barcodes)

combined_metadata <- all_except_subsets_metadata |>
  bind_rows(T_V2_metadata) |>
  bind_rows(B_V2_metadata) |>
  # filter(!(barcode %in% T_contaminating_barcodes))
nrow(combined_metadata)
combined_metadata |> head(3)

combined_metadata |>
  write_tsv(file.path(out_dir, 'combined_metadata.tsv'))
  # write_tsv(file.path(out_dir, 'combined_metadata_cleaned.tsv'))
