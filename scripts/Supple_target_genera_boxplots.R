# ==============================================================================
# Script: Supple_target_genera_boxplots.R
# Author: Kotaro Murai
# Description: 
#   Extract specific genera (Prevotella, Succinivibrio, Treponema) from mOTUs3 data
#   by aggregating all sub-lineages (e.g., Treponema_A, Treponema_G).
# ==============================================================================

# Load libraries
library(tidyverse)
library(ggpubr)
library(patchwork)

# Load custom plot styles
source("scripts/00_plot_styles.R")

# Create output directories
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Data Loading 
# ==============================================================================
cat("Loading data...\n")
Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

# Load metadata
metadata <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)

# Enforce the group order: Africa -> Non-Ind -> Ind
metadata$Analysis_Group <- factor(
  metadata$Analysis_Group, 
  levels = c("Sub-Saharan Africa", "Other Non-Industrialized", "Industrialized")
)

# Load abundance table (mOTUs3)
species_data <- read_tsv("data/metalog_motus3_processed.tsv.gz", show_col_types = FALSE) %>%
  column_to_rownames(var = "rowname")

# Load taxonomy
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"

# Create a mapping table between column names and taxonomy
column_mapping <- data.frame(original_col = colnames(species_data)) %>%
  mutate(id = str_extract(original_col, "(?<=\\[)[^\\[\\]]+(?=\\]$)")) %>%
  left_join(taxa_data, by = "id")

# ==============================================================================
# 2. Extract & Aggregate Target Genera
# ==============================================================================
targets <- c("Prevotella", "Succinivibrio", "Treponema")
target_abundances <- list()

cat("Aggregating targeted genera...\n")

for (target_genus in targets) {
  # Find all columns where the genus string contains the target name
  target_motus <- column_mapping %>% 
    filter(str_detect(genus, target_genus)) %>% 
    pull(original_col)
  
  if(length(target_motus) > 0) {
    if(length(target_motus) == 1) {
      aggregated_abun <- species_data[[target_motus]]
    } else {
      aggregated_abun <- rowSums(species_data[, target_motus, drop = FALSE], na.rm = TRUE)
    }
  } else {
    aggregated_abun <- rep(0, nrow(species_data))
    cat(sprintf("Warning: No mOTUs found for %s\n", target_genus))
  }
  
  target_abundances[[target_genus]] <- aggregated_abun
}

# Combine into a single dataframe and merge with metadata
df_aggregated <- bind_cols(target_abundances)
rownames(df_aggregated) <- rownames(species_data)

plot_data <- df_aggregated %>%
  rownames_to_column("sample_alias") %>%
  left_join(metadata, by = "sample_alias") %>%
  filter(!is.na(Analysis_Group)) %>%
  pivot_longer(cols = all_of(targets), names_to = "Genus", values_to = "abundance") %>%
  mutate(
    # Log10 transform to match the assumption of the LMM analysis
    log10_abundance = log10(abundance + 1E-4)
  )

# ==============================================================================
# 3. Plotting
# ==============================================================================
cat("Generating plots...\n")

# Define comparison pairs based on the factor levels
my_comparisons <- list(
  c("Sub-Saharan Africa", "Other Non-Industrialized"),
  c("Other Non-Industrialized", "Industrialized"),
  c("Sub-Saharan Africa", "Industrialized")
)

global_ymin <- min(plot_data$log10_abundance, na.rm = TRUE)
global_ymax <- max(plot_data$log10_abundance, na.rm = TRUE)
step_size <- (global_ymax - global_ymin) * 0.1
y_pos <- c(
  global_ymax + step_size * 1, 
  global_ymax + step_size * 2,
  global_ymax + step_size * 3 
)
plot_ymax <- y_pos[3] + step_size
# Custom function to create polished boxplots for each genus
create_genus_boxplot <- function(data, genus_name) {
  
  g_data <- data %>% filter(Genus == genus_name)
  if (genus_name == "Treponema") {
    y_label <- bquote(log[10]~italic(.(genus_name)~"sensu lato")~"abundance")
  } else {
    y_label <- bquote(log[10]~italic(.(genus_name))~"abundance")
  }
  
  ggplot(g_data, aes(x = Analysis_Group, y = log10_abundance, fill = Analysis_Group)) +
    geom_jitter(width = 0.2, alpha = 0.1, size = 0.05, color = "grey50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black", lwd = 0.25) +
    
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      p.adjust.method = "bonferroni",
      label = "p.signif",
      vjust = 0.5,
      size = 2,                      
      tip.length = 0.01,
      y.position = y_pos
    ) +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    scale_y_continuous(
      limits = c(global_ymin, plot_ymax),
      expand = expansion(mult = c(0.05, 0))
    ) +
    labs(x = NULL, y = y_label) +
    theme_nature(base_size = 6) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 5.5, lineheight = 0.9), 
      panel.grid.major.x = element_blank()
    )
}

# Generate individual plots
plots_list <- map(targets, ~create_genus_boxplot(plot_data, .x))

# Combine plots vertically using patchwork
p_combined <- wrap_plots(plots_list, nrow = 1)

# ==============================================================================
# 4. Save Output
# ==============================================================================
filename_base <- "results/figures/Supplementary_Figure1_Target_Genera_Boxplots"

ggsave(paste0(filename_base, ".pdf"), p_combined, width = 180, height = 60, units = "mm", useDingbats = FALSE)
ggsave(paste0(filename_base, ".png"), p_combined, width = 180, height = 60, units = "mm", bg = "white", dpi = 300)

cat(sprintf("Successfully saved plots to %s(.pdf/.png)\n", filename_base))