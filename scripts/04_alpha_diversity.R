# ==============================================================================
# Script: 04_alpha_diversity.R
# Author: Kotaro Murai
# Description: 
#   Calculate Alpha Diversity (Shannon index, Observed species).
#   Generate boxplots with Wilcoxon pairwise comparisons (Bonferroni adjusted).
# ==============================================================================

# Load libraries
library(tidyverse)
library(vegan)
library(ggpubr)
library(patchwork)

# Load custom plot styles
source("scripts/00_plot_styles.R")

# Create output directories if they don't exist
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---
cat("Loading data...\n")
Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

# Load abundance table and set row names to sample IDs for vegan
mat <- read_tsv("data/metalog_motus3_processed.tsv.gz") %>%
  column_to_rownames(var = "rowname")

# Load pre-grouped metadata from Step 02
metadata <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)

# Enforce the group order: Africa -> Non-Ind -> Ind
metadata$Analysis_Group <- factor(
  metadata$Analysis_Group, 
  levels = c("Africa", "Non-Africa Non-Industrialized", "Industrialized")
)

# --- 2. Calculate Alpha Diversity ---
cat("Calculating Alpha Diversity...\n")

# Calculate diversity indices
alpha_div <- tibble(
  sample_alias = rownames(mat),
  shannon = diversity(mat, index = "shannon"), # Shannon Index
  observed = specnumber(mat)                   # Richness (Observed species)
)

# Merge with metadata
# Note: Ensure the joining key matches. Assuming 'sample_alias' is the key.
df_analysis <- metadata %>%
  left_join(alpha_div, by = "sample_alias")

# --- 3. Plotting Setup ---
cat("Plotting...\n")

# Define comparison pairs based on the factor levels
my_comparisons <- list(
  c("Africa", "Non-Africa Non-Industrialized"),
  c("Non-Africa Non-Industrialized", "Industrialized"),
  c("Africa", "Industrialized")
)

# Custom function to create polished boxplots
create_boxplot <- function(data, y_col, y_label, title_text) {
  ggplot(data, aes(x = Analysis_Group, y = .data[[y_col]], fill = Analysis_Group)) +
    
    # Boxplot
    geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black", lwd = 0.4) +
    
    # Jittered points for individual samples
    geom_jitter(width = 0.2, alpha = 0.1, size = 0.2) +
    
    # Statistical testing: Wilcoxon rank-sum test with Bonferroni correction
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      p.adjust.method = "bonferroni", # 論文テキストとの整合性のために追加
      label = "p.signif",             # Asterisks
      vjust = 0.5,
      size = 5,                       
      tip.length = 0.01
    ) +
    
    # Apply global colors from 00_plot_styles.R
    scale_fill_manual(values = thesis_colors) +
    
    # Wrap x-axis labels
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) +
    
    # Expand upper y-axis margin to prevent asterisk cutoff
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    
    labs(x = NULL, y = y_label, title = title_text) +
    
    # Apply global theme and override specific elements
    theme_thesis(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(size = 11, lineheight = 0.9), 
      axis.text.y = element_text(size = 11),
      panel.grid.major.x = element_blank()
    )
}

# Generate individual plots
p1 <- create_boxplot(df_analysis, "observed", "Observed Species", "Richness (Observed Species)")
p2 <- create_boxplot(df_analysis, "shannon", "Shannon Index", "Shannon Index")

# Combine plots using patchwork
p_combined <- p1 + p2

# --- 4. Save Output ---
filename_base <- "results/figures/Figure1d_Alpha_Diversity"

ggsave(paste0(filename_base, ".pdf"), p_combined, width = 10, height = 6, useDingbats = FALSE)
ggsave(paste0(filename_base, ".png"), p_combined, width = 10, height = 6, dpi = 300, bg = "white")

print(p_combined)
cat(sprintf("Successfully saved plots to %s(.pdf/.png)\n", filename_base))