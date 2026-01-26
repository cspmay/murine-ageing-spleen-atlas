library(VennDiagram)
library(ggplot2)
library(dplyr)
library(readr)
library(grid)
library(gridExtra)

out_dir <- './NM_mice/out/pseudobulk_DE'

### TWO-WAY VENN ###
create_venn_diagram_two_cell_types <- function(cell_type_1, cell_type_2, de_results_dir, output_dir = NULL) {
  file_1 <- file.path(de_results_dir, paste0(cell_type_1, '_age_sig.tsv'))
  file_2 <- file.path(de_results_dir, paste0(cell_type_2, '_age_sig.tsv'))  
  sig_1 <- read_tsv(file_1, show_col_types = FALSE)
  sig_2 <- read_tsv(file_2, show_col_types = FALSE) 
  all_sig_genes_1 <- sig_1 %>% pull(name)
  all_sig_genes_2 <- sig_2 %>% pull(name) 
  shared_genes <- inner_join(sig_1, sig_2, by = "name", suffix = c("_ct1", "_ct2"))
  shared_genes <- shared_genes %>%
    mutate(
      same_direction = sign(lfc_ct1) == sign(lfc_ct2),
      regulation = case_when(
        lfc_ct1 > 0 & lfc_ct2 > 0 ~ "upregulated",
        lfc_ct1 < 0 & lfc_ct2 < 0 ~ "downregulated",
        TRUE ~ "discordant"
      )
    )  
  shared_genes_concordant <- shared_genes %>% filter(same_direction)
  upreg_genes_1 <- sig_1 %>% filter(lfc > 0.25) %>% pull(name)
  upreg_genes_2 <- sig_2 %>% filter(lfc > 0.25) %>% pull(name)
  downreg_genes_1 <- sig_1 %>% filter(lfc < -0.25) %>% pull(name)
  downreg_genes_2 <- sig_2 %>% filter(lfc < -0.25) %>% pull(name) 
  common_upreg <- shared_genes_concordant %>% 
    filter(regulation == "upregulated") %>% 
    pull(name)
  common_downreg <- shared_genes_concordant %>% 
    filter(regulation == "downregulated") %>% 
    pull(name)
  if (is.null(output_dir)) {
    shared_genes_path <- file.path(de_results_dir, 'shared_genes')
  } else {
    shared_genes_path <- output_dir
  }
  dir.create(shared_genes_path, showWarnings = FALSE)
  venn_all <- venn.diagram(
    x = setNames(
    list(all_sig_genes_1, all_sig_genes_2),
    c(cell_type_1, cell_type_2)
    ),
    category.names = c(
      paste0(cell_type_1, "\n(", length(all_sig_genes_1), ")"),
      paste0(cell_type_2, "\n(", length(all_sig_genes_2), ")")
    ),
    filename = NULL,
    fill = c("#ff6b6b", "#4ecdc4"),
    alpha = 0.6,
    cex = 2.0, 
    cat.cex = 2.0, 
    cat.pos = c(30, -30),
    euler.d = FALSE,  
    scaled = FALSE
  )
  venn_filename <- file.path(shared_genes_path, 
                             paste0("venn_diagram_", cell_type_1, "_vs_", cell_type_2, ".png"))
  png(
    filename = venn_filename,
    width = 8,
    height = 8,
    units = "in",
    res = 600
  )
  grid.draw(venn_all)
  dev.off()
}

### THREE-WAY VENN ###
create_venn_diagram_three_cell_types <- function(cell_type_1, cell_type_2, cell_type_3, de_results_dir, output_dir = NULL) {
  file_1 <- file.path(de_results_dir, paste0(cell_type_1, '_age_sig.tsv'))
  file_2 <- file.path(de_results_dir, paste0(cell_type_2, '_age_sig.tsv'))
  file_3 <- file.path(de_results_dir, paste0(cell_type_3, '_age_sig.tsv'))
  sig_1 <- read_tsv(file_1, show_col_types = FALSE)
  sig_2 <- read_tsv(file_2, show_col_types = FALSE)
  sig_3 <- read_tsv(file_3, show_col_types = FALSE)
  all_sig_genes_1 <- sig_1 %>% pull(name)
  all_sig_genes_2 <- sig_2 %>% pull(name)
  all_sig_genes_3 <- sig_3 %>% pull(name)
  upreg_genes_1 <- sig_1 %>% filter(lfc > 0.25) %>% pull(name)
  upreg_genes_2 <- sig_2 %>% filter(lfc > 0.25) %>% pull(name)
  upreg_genes_3 <- sig_3 %>% filter(lfc > 0.25) %>% pull(name)
  downreg_genes_1 <- sig_1 %>% filter(lfc < -0.25) %>% pull(name)
  downreg_genes_2 <- sig_2 %>% filter(lfc < -0.25) %>% pull(name)
  downreg_genes_3 <- sig_3 %>% filter(lfc < -0.25) %>% pull(name)
  common_upreg_all <- Reduce(intersect, list(upreg_genes_1, upreg_genes_2, upreg_genes_3))
  common_upreg_1_2 <- setdiff(intersect(upreg_genes_1, upreg_genes_2), upreg_genes_3)
  common_upreg_1_3 <- setdiff(intersect(upreg_genes_1, upreg_genes_3), upreg_genes_2)
  common_upreg_2_3 <- setdiff(intersect(upreg_genes_2, upreg_genes_3), upreg_genes_1)
  common_downreg_all <- Reduce(intersect, list(downreg_genes_1, downreg_genes_2, downreg_genes_3))
  common_downreg_1_2 <- setdiff(intersect(downreg_genes_1, downreg_genes_2), downreg_genes_3)
  common_downreg_1_3 <- setdiff(intersect(downreg_genes_1, downreg_genes_3), downreg_genes_2)
  common_downreg_2_3 <- setdiff(intersect(downreg_genes_2, downreg_genes_3), downreg_genes_1)
  if (is.null(output_dir)) {
    shared_genes_path <- file.path(de_results_dir, 'shared_genes')
  } else {
    shared_genes_path <- output_dir
  }
  dir.create(shared_genes_path, showWarnings = FALSE)
  venn_all <- venn.diagram(
    x = list(
      all_sig_genes_1,
      all_sig_genes_2,
      all_sig_genes_3
    ),
    category.names = c(
      paste0(cell_type_1, "\n(", length(all_sig_genes_1), ")"),
      paste0(cell_type_2, "\n(", length(all_sig_genes_2), ")"),
      paste0(cell_type_3, "\n(", length(all_sig_genes_3), ")")
    ),
    filename = NULL,
    fill = c("#ff6b6b", "#4ecdc4", "#95e1d3"),
    alpha = 0.5,
    cex = 1.5,  
    cat.cex = 1.5,  
    cat.pos = c(-27, 27, 180),  
    cat.dist = c(0.075, 0.075, 0.065)
  )
  venn_filename <- file.path(shared_genes_path, 
                             paste0("venn_diagram_", cell_type_1, "_", cell_type_2, "_", cell_type_3, ".png"))
  png(
    filename = venn_filename,
    width = 8,
    height = 8,
    units = "in",
    res = 600
  ) 
  grid.draw(venn_all)
  dev.off()
}


### T CELLS ###
deseq2_out_dir <- file.path(out_dir, 'T', 'age')

cell_type_1 <- "CD4N"
cell_type_2 <- "CD8N"
overlap_genes <- create_venn_diagram_two_cell_types(cell_type_1, cell_type_2, deseq2_out_dir)

cell_type_1 <- "CD4EM"
cell_type_2 <- "CD8EM"
cell_type_3 <- "CD8EM Gzmk+"
overlap_genes_3way <- create_venn_diagram_three_cell_types(cell_type_1, cell_type_2, cell_type_3, deseq2_out_dir)


### B CELLS ###
cell_type_1 <- "MZB"
cell_type_2 <- "TrB"
cell_type_3 <- "FoB"
overlap_genes_3way <- create_venn_diagram_three_cell_types(cell_type_1, cell_type_2, cell_type_3, deseq2_out_dir)

cell_type_1 <- "Memory B"
cell_type_2 <- "Apoe-hi B"
cell_type_3 <- "ABCs"
overlap_genes_3way <- create_venn_diagram_three_cell_types(cell_type_1, cell_type_2, cell_type_3, deseq2_out_dir)
