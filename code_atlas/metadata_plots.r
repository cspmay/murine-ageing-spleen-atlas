library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

out_dir <- './atlas/out/'
out_dir

plot_dir = './out/figures'
plot_dir

combined_metadata <- file.path(out_dir, 'combined_metadata.tsv') |>
    read_delim(show_col_types = F)  |>
    select(barcode, sample_name, age_group, condition, dataset, cell_type_V2)
combined_metadata |> head(3)
unique(combined_metadata$cell_type_V2)

### CELL DISTRIBUTION ACROSS CONDITIONS ###
shorten_number <- function(x) {
  if (x >= 1e6) {
    return(paste0(round(x / 1e6, 1), "M"))
  } else if (x >= 1e3) {
    return(paste0(round(x / 1e3, 1), "k"))
  } else {
    return(as.character(x))
  }
}
combined_metadata <- combined_metadata |> 
  mutate(
    condition = recode(condition, 
      "Het" = "Heterochronic Parabiosis",
      "Iso" = "Isochronic Parabiosis",
      "PF4" = "PF4 Treatment",
      "control" = "Control",
      "stress" = "Stress"
    ),
    age_group = recode(age_group,
      "old" = "Old",
      "young" = "Young"
    )
  )
p <- combined_metadata |> 
  count(dataset, age_group, condition, name = "num_cells") |>
  mutate(label = sapply(num_cells, shorten_number)) |>  
  ggplot(aes(x = num_cells, y = dataset, fill = interaction(age_group, condition))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4, color = "black") + 
  scale_fill_brewer(palette = "Set2") +  
  theme_minimal() +
  labs(
    x = "Number of Cells",
    y = "Dataset",
    fill = "Age Group & Condition"
  )
ggsave(
  file = file.path(plot_dir, 'all_cells.png'),
  plot = p,
  device = 'png',
  dpi = 600
)
p

### CELL DISTRIBUTION ACROSS AGE ###
plot_data <- combined_metadata |> 
  count(dataset, age_group, condition, name = "num_cells") |>
  filter(condition %in% c('Control', 'Isochronic Parabiosis')) |>
  group_by(dataset) |>
  mutate(
    total_cells = sum(num_cells),
    label = sapply(total_cells, shorten_number)
  ) |>
  ungroup()
p <- plot_data |>
  ggplot(aes(x = num_cells, y = dataset, fill = age_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(x = total_cells, label = label), 
            hjust = -0.2, size = 5, color = "black") +  
  scale_x_continuous(limits = c(0, 60000), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2") +  
  theme_minimal() +
  labs(
    x = "Number of Cells",
    y = "Dataset",
    fill = "Age Group"
  )
ggsave(
  file = file.path(plot_dir, 'ageing_Y_O_cells.png'),
  plot = p,
  device = 'png',
  width = 14,
  height = 6,
  dpi = 600
)
p