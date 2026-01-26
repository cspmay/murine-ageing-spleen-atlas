library(tidyverse)
library(scran)
library(Matrix)
library(glmGamPoi)
library(BiocSingular)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(DESeq2)

init_dir <- function(new_dir) {
  if (!dir.exists(new_dir)) {
	dir.create(new_dir)
  }
}

out_dir <- './NM_mice/out/pseudobulk_DE'

DE_out_dir <- file.path(out_dir, 'all')
# DE_out_dir <- file.path(out_dir, 'T')
# DE_out_dir <- file.path(out_dir, 'B')
# DE_out_dir <- file.path(out_dir, 'Gzmk')
init_dir(DE_out_dir)
print(DE_out_dir)

load_sce <- function(sce_path) {
  sce <- sce_path |>
    readRDS() |>
    prepare_sce()  
  return(sce)
}

prepare_sce <- function(sce) {
  result_sce <- logNormCounts(sce)
  top_features <- getTopHVGs(result_sce, n = 5000)
  result_sce <- fixedPCA(
    result_sce,
    subset.row = top_features, 
    BSPARAM=RandomParam()
  )
  return(result_sce)
}

NM_mice_sce <- file.path(out_dir, 'NM_mice_all_stress_pseudobulk.rds') |>
# NM_mice_sce <- file.path(out_dir, 'NM_mice_T_stress_pseudobulk.rds') |>
# NM_mice_sce <- file.path(out_dir, 'NM_mice_B_stress_pseudobulk.rds') |>
# NM_mice_sce <- file.path(out_dir, 'NM_mice_Gzmk_stress_pseudobulk.rds') |>
  load_sce()

pca_coords <- reducedDim(NM_mice_sce, 'PCA')
pca_data <- data.frame(
  PC1 = pca_coords[, 1],
  PC2 = pca_coords[, 2],
  condition = colData(NM_mice_sce)$condition,
  cell_type_V2 = colData(NM_mice_sce)$cell_type_V2
)
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(alpha = 0.6, size = 0.5) +
  theme_minimal() +
  labs(title = "PCA - Condition")
p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_type_V2)) +
  geom_point(alpha = 0.6, size = 0.5) +
  theme_minimal() +
  labs(title = "PCA - Cell Type") +
  theme(legend.position = "right")
p1 + p2

### DESeq2 ###
run_DE_deseq2 <- function(sce, out_dir, batch_name) {
  batch_out_dir <- file.path(out_dir, batch_name)
  init_dir(batch_out_dir) 
  batch_plots_out_dir <- file.path(batch_out_dir, 'plots')
  init_dir(batch_plots_out_dir)
  count_matrix <- counts(sce)
  col_data <- colData(sce) |> as.data.frame()
  col_data$condition <- factor(col_data$condition, levels = c("control", "stress"))  
  for(current_cell_type in unique(col_data$cell_type_V2)) {
    cell_type_idx <- col_data$cell_type_V2 == current_cell_type
    cell_type_counts <- count_matrix[, cell_type_idx]
    cell_type_coldata <- col_data[cell_type_idx, ]
    n_control <- sum(cell_type_coldata$condition == "control")
    n_stress <- sum(cell_type_coldata$condition == "stress")
    dds <- DESeqDataSetFromMatrix(
      countData = round(cell_type_counts),
      colData = cell_type_coldata,
      design = ~ condition
    )
    dds <- dds[rowSums(counts(dds)) >= 10, ]   
    dds <- DESeq(dds, quiet = TRUE)   
    res <- results(dds, contrast = c("condition", "stress", "control"))    
    res_df <- as.data.frame(res) |>
      tibble::rownames_to_column("name") |>
      dplyr::rename(
        lfc = log2FoldChange,
        pval = pvalue,
        adj_pval = padj
      ) |>
      dplyr::filter(!is.na(adj_pval)) |>
      dplyr::arrange(lfc)
    cell_type_for_file <- gsub('/', '_', current_cell_type)
    res_df |>
      dplyr::filter(adj_pval < 0.05, abs(lfc) > 0.25) |>
      write_tsv(file.path(batch_out_dir, paste0(cell_type_for_file, '_stress_sig.tsv')))
  }
}

run_DE_deseq2(NM_mice_sce, DE_out_dir, 'stress')


### VOLCANO PLOTS ###
deseq2_out_dir <- file.path(DE_out_dir, 'stress')

create_volcano_plots <- function(deseq2_out_dir, plot_dir) {
  plot_dir <- file.path(deseq2_out_dir, 'publication_plots')
  init_dir(plot_dir)
  result_files <- list.files(
    deseq2_out_dir,
    pattern = "_stress_all.tsv$",
    full.names = TRUE
  )  
  for (file in result_files) {
    filename <- basename(file)
    cell_type_name <- sub("_stress_all.tsv$", "", filename)
    cell_type <- gsub("_", " ", cell_type_name)
    result <- read_tsv(file, show_col_types = FALSE)
    res_sig <- result |>
      filter(adj_pval < 0.05, abs(lfc) > 0.25)
    n_up <- sum(res_sig$lfc > 0)
    n_down <- sum(res_sig$lfc < 0)
    plot_title <- paste0(cell_type)  
    plot_data <- result |>
      filter(!is.na(adj_pval)) |>
      mutate(
        is_de = case_when(
          adj_pval < 0.05 & lfc > 0.25 ~ 'Upregulated in stress',
          adj_pval < 0.05 & lfc < -0.25 ~ 'Downregulated in stress',
          TRUE ~ 'Not changed'
        ),
        label = ifelse(adj_pval < 0.05 & abs(lfc) > 0.25, name, NA)
      )    
    color_values <- c(
      'Downregulated in stress' = '#415484', 
      'Upregulated in stress' = '#D55640',
      'Not changed' = 'grey60'
    )
    p <- ggplot(plot_data, aes(x = lfc, y = -log10(adj_pval), label = label, col = is_de)) +
      geom_point(size = 2, alpha = 0.8) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = 'dashed', alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = 'dashed', alpha = 0.5) +
      geom_label_repel(
        show.legend = FALSE, 
        box.padding = 0.25, 
        size = 4, 
        max.overlaps = 25, 
        min.segment.length = Inf
      ) +
      scale_color_manual(values = color_values) +
      theme_minimal() +
      theme(
        legend.position = 'top', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
      ) +
      guides(col = guide_legend(title = NULL)) +
      labs(
        x = expression(log[2]~"(fold change)"), 
        y = expression(-log[10]~"(adjusted"~italic(p)*"-value)"),
        title = plot_title
      )  
    ggsave(
      file = file.path(plot_dir, paste0(cell_type_name, "_stress.png")),
      plot = p,
      device = 'png',
      dpi = 600
    )

  }
}

create_volcano_plots(deseq2_out_dir, plot_dir)
