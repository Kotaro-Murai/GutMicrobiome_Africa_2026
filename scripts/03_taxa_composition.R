# ==============================================================================
# Script: 03_taxa_composition.R
# Author: Kotaro Murai
# Description: 
#   Generate taxonomic composition barplots (Figure 1c) at a specified rank.
#   Samples are aggregated by country and separated by Analysis Group.
# ==============================================================================

# Load libraries
library(tidyverse)
library(data.table)
library(patchwork)

# Load custom plot styles and palettes (including 'taxa_colors_20')
source("scripts/00_plot_styles.R")

# Create output directories if they don't exist
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---
cat("Loading data...\n")
Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

# Load species abundance data
species_data <- read_tsv("data/metalog_motus3_processed.tsv.gz") %>%
  column_to_rownames(var = "rowname")
# Load taxonomic annotation
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"

# Load the metadata that already contains the 'Analysis_Group' assigned in Step 02
metadata <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)

metadata$Analysis_Group <- factor(metadata$Analysis_Group, 
                                  levels = c("Africa", "Non-Africa Non-Industrialized", "Industrialized"))

# Remove noise column if exists
if("-1" %in% colnames(species_data)) {
  species_data <- species_data %>% select(-`-1`)
}

# --- 2. Analysis Settings ---
target_rank <- "phylum"  # Can be changed to "genus", "order", etc.
top_n_taxa <- 10         # Number of top taxa to display

# ==============================================================================
# 3. Data Aggregation and Formatting
# ==============================================================================
cat("Aggregating taxa data at", target_rank, "level...\n")

# Create mapping between column names and GTDB taxonomy
col_names <- colnames(species_data)
ids <- str_extract(col_names, "(?<=\\[)[^\\[\\]]+(?=\\]$)")

map_df <- data.frame(original_col = col_names, id = ids) %>%
  left_join(taxa_data, by = "id") %>%
  mutate(
    RankName = .[[target_rank]],
    RankName = if_else(is.na(RankName) | RankName == "", "Unassigned", RankName),
    CleanName = str_remove(RankName, "^[pcofg]__") # Remove GTDB prefixes
  ) %>%
  filter(!is.na(id))

# Aggregate abundance by the specified taxonomic rank
abd_mat <- as.matrix(species_data[, map_df$original_col])
unique_taxa <- unique(map_df$CleanName)

rank_abd_list <- lapply(unique_taxa, function(tax) {
  cols <- map_df$original_col[map_df$CleanName == tax]
  if(length(cols) == 1) return(abd_mat[, cols])
  return(rowSums(abd_mat[, cols]))
})
rank_abd_df <- as.data.frame(do.call(cbind, rank_abd_list))
colnames(rank_abd_df) <- unique_taxa

# Calculate relative abundance per country, keeping Analysis_Group
plot_source <- bind_cols(metadata %>% select(country, Analysis_Group), rank_abd_df) %>%
  pivot_longer(cols = -c(country, Analysis_Group), names_to = "Taxon", values_to = "Abundance") %>%
  group_by(country, Analysis_Group, Taxon) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(country) %>%
  mutate(RelAbundance = MeanAbundance / sum(MeanAbundance)) %>%
  ungroup()

# --- 4. Handle "Others" Category ---
# Identify the top N most abundant taxa overall
top_taxa_list <- plot_source %>%
  filter(Taxon != "Unassigned") %>%
  group_by(Taxon) %>%
  summarise(TotalMean = mean(RelAbundance)) %>%
  arrange(desc(TotalMean)) %>%
  slice_head(n = top_n_taxa) %>%
  pull(Taxon)

# Define factor levels (Bacteroides/Prevotella first if present, then top taxa, then Others)
custom_levels <- unique(c("Bacteroides", "Prevotella", top_taxa_list))
custom_levels <- custom_levels[custom_levels %in% top_taxa_list]
custom_levels <- c(custom_levels, "Others") 

plot_data <- plot_source %>%
  mutate(
    DisplayTaxon = if_else(Taxon %in% top_taxa_list, Taxon, "Others"),
    DisplayTaxon = factor(DisplayTaxon, levels = custom_levels)
  ) %>%
  group_by(country, Analysis_Group, DisplayTaxon) %>%
  summarise(RelAbundance = sum(RelAbundance), .groups = "drop")

# ==============================================================================
# 5. Sorting Control (Sort countries by the dominant taxon in their group)
# ==============================================================================
cat("Calculating sort order per group...\n")

# 1. Identify the most abundant taxon within each Analysis Group
group_dominant_taxa <- plot_data %>%
  filter(DisplayTaxon != "Others") %>%
  group_by(Analysis_Group, DisplayTaxon) %>%
  summarise(MeanAb = mean(RelAbundance), .groups = "drop") %>%
  group_by(Analysis_Group) %>%
  slice_max(MeanAb, n = 1) %>% 
  select(Analysis_Group, DominantTaxon = DisplayTaxon)

print(group_dominant_taxa)

# 2. Sort countries based on the abundance of their group's dominant taxon
country_order <- plot_data %>%
  left_join(group_dominant_taxa, by = "Analysis_Group") %>%
  filter(DisplayTaxon == DominantTaxon) %>% 
  arrange(Analysis_Group, desc(RelAbundance)) %>% 
  pull(country)

# 3. Apply factor levels to country names to enforce sorting
plot_data$country <- factor(plot_data$country, levels = country_order)

# --- 6. Plotting ---
cat("Plotting...\n")

plot_taxa <- levels(plot_data$DisplayTaxon)
my_colors <- character(length(plot_taxa))
names(my_colors) <- plot_taxa

fallback_idx <- 1 

for (tax in plot_taxa) {
  if (tax %in% names(phylum_master_colors)) {
    my_colors[tax] <- phylum_master_colors[tax]
  } else if (tax == "Others") {
    my_colors[tax] <- "grey90"
  } else {
    my_colors[tax] <- taxa_colors_20[fallback_idx]
    fallback_idx <- fallback_idx + 1
  }
}

p_bar <- ggplot(plot_data, aes(x = country, y = RelAbundance, fill = DisplayTaxon)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.1) + 
  
  scale_fill_manual(values = my_colors, name = NULL) + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  
  # Split by Analysis Group
  facet_grid(~ Analysis_Group, scales = "free_x", space = "free_x") +
  
  labs(
    title = paste0("Gut Microbiota Composition (Rank: ", str_to_title(target_rank), ")"),
    y = "Relative Abundance",
    x = NULL
  ) +
  
  # Apply base theme but override specific elements for this barplot
  theme_thesis(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    strip.text = element_text(face = "bold", size = 14, color = "white"),
    strip.background = element_rect(fill = "#555555", color = NA) 
  )

# --- 7. Save Output ---
filename_base <- paste0("results/figures/Figure1c_Composition_", target_rank)

# Save as PDF for Illustrator (Vector) and PNG for quick preview
ggsave(paste0(filename_base, ".pdf"), p_bar, width = 18, height = 9, useDingbats = FALSE)
ggsave(paste0(filename_base, ".png"), p_bar, width = 18, height = 9, dpi = 300, bg = "white")

print(p_bar)
cat(sprintf("Successfully saved plots to %s(.pdf/.png)\n", filename_base))