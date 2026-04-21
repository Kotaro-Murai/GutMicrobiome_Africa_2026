# ==============================================================================
# Script: 06_Lollipop_Species.R
# Description: Species-level Lollipop Plot (Figure 2b).
# ==============================================================================

library(tidyverse)
library(data.table)

# Load custom plot styles
source("scripts/00_plot_styles.R")

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load Data
# ==============================================================================
cat("Loading LMM species results and taxonomy...\n")

# Load the wide format results from Script 05
results_wide <- read_csv("results/tables/lmm_results_species_wide.csv", show_col_types = FALSE)

# Load taxonomy to get Phylum information (for consistent coloring with Fig 2c)
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"

# ==============================================================================
# 2. Data Preparation (Dual-filtering & Conservative Estimate)
# ==============================================================================
cat("Preparing data for plotting...\n")

# Extract mOTU ID and merge with taxonomy
plot_data <- results_wide %>%
  mutate(mOTU_id = str_extract(species, "(?<=\\[)[^\\[\\]]+(?=\\][^\\[]*$)")) %>%
  left_join(taxa_data, by = c("mOTU_id" = "id")) %>%
  filter(sign(estimate_SSAvsInd) == sign(estimate_SSAvsONI)) %>%
  filter(q.value_SSAvsInd < 0.05 & q.value_SSAvsONI < 0.05) %>%
  mutate(
    conservative_est = if_else(
      abs(estimate_SSAvsInd) < abs(estimate_SSAvsONI),
      estimate_SSAvsInd,
      estimate_SSAvsONI
    ),
    estimate = conservative_est,
    Type = if_else(estimate > 0, "Enriched", "Depleted"),
    
    # 4. Clean up names: Remove [mOTU...], s__, and Unassigned
    clean_name = coalesce(species.y, species.x) %>% 
      str_remove("\\[[^\\[\\]]+\\]$") %>% 
      str_remove("^s__") %>% 
      str_trim(),
    clean_phylum = str_remove(phylum, "^p__")
  ) %>%
  filter(!str_detect(clean_name, "Not_annotated|Incongruent|-1|Unassigned"))

plot_data %>%
  arrange(desc(estimate)) %>%
  write_csv("results/tables/Supplementary_Table_SSA_Specific_Species.csv")

# Selection Logic: Top 10 Enriched and Top 10 Depleted
top_n_show <- 10
est_threshold <- 0.1

plot_enriched <- plot_data %>%
  filter(Type == "Enriched", estimate > est_threshold) %>%
  arrange(desc(estimate)) %>% slice_head(n = top_n_show)

plot_depleted <- plot_data %>%
  filter(Type == "Depleted", estimate < -est_threshold) %>%
  arrange(estimate) %>% slice_head(n = top_n_show)

final_plot_data <- bind_rows(plot_enriched, plot_depleted) %>%
  arrange(estimate) %>%
  distinct(species.x, .keep_all = TRUE) %>%
  mutate(
    # Ensure factor levels follow the estimate for correct Y-axis ordering
    unique_key = factor(species.x, levels = unique(species.x))
  )

# ==============================================================================
# 3. Apply Master Colors (Consistent with Script 07)
# ==============================================================================
# Status colors (Africa vs Non-Africa)
cols_status <- c("Enriched" =LMM_colors[["Enriched"]], 
                 "Depleted" = LMM_colors[["Depleted"]])

# Phylum colors
used_phyla <- unique(final_plot_data$clean_phylum)
defined_cols <- phylum_master_colors[names(phylum_master_colors) %in% used_phyla]
missing_phyla <- setdiff(used_phyla, names(phylum_master_colors))

if(length(missing_phyla) > 0) {
  extra_cols <- taxa_colors_20[1:length(missing_phyla)]
  names(extra_cols) <- missing_phyla
  final_phylum_colors <- c(defined_cols, extra_cols)
} else {
  final_phylum_colors <- defined_cols
}

# ==============================================================================
# 4. Visualization
# ==============================================================================
cat("Generating Figure 2b (Species Lollipop)...\n")

y_labels <- setNames(final_plot_data$clean_name, final_plot_data$unique_key)

p <- ggplot(final_plot_data, aes(x = estimate, y = unique_key)) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.25) +
  geom_vline(xintercept = c(est_threshold, -est_threshold), linetype = "dashed", color = "gray80", linewidth = 0.25) +
  
  # Stick (Segment)
  geom_segment(aes(x = 0, xend = estimate, y = unique_key, yend = unique_key, color = Type), linewidth = 0.5) +
  scale_color_manual(values = cols_status, guide = "none") + 
  
  # Head (Point) - Using Shape 21 (Circle with border) to show Phylum fill
  geom_point(aes(fill = clean_phylum), shape = 21, size = 2, color = "black", stroke = 0.25) + 
  scale_fill_manual(values = final_phylum_colors, name = "Phylum", guide = "none") +
  
  # Formatting
  scale_y_discrete(labels = y_labels) +
  scale_x_continuous(expand = expansion(mult = 0.15)) + 
  labs(x = "LMM Estimate", y = NULL, title = NULL) +
  
  theme_nature(base_size = 6) +
  theme(
    axis.text.y = element_text(face = "italic", size = 7, color = "black"), 
    legend.position = "bottom", 
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.25)
  )

# Save
outfile <- "results/figures/Figure2b_Lollipop_Species"
ggsave(paste0(outfile, ".pdf"), p, width = 88, height = 80, units = "mm", useDingbats = FALSE, bg = "transparent")
ggsave(paste0(outfile, ".png"), p, width = 88, height = 80, units = "mm", dpi = 300, bg = "white")

cat("Successfully saved Figure 2b!\n")