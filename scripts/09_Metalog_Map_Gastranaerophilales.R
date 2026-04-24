# ==============================================================================
# Script: 09_Metalog_map_Gastranaerophilales.R
# Description: Generates map and bar plot (Figure 3a) of 
#              Gastranaerophilales prevalence using Metalog (Metagenome) data.
# Author: Kotaro Murai
# ==============================================================================

library(tidyverse)
library(data.table)
library(sf)
library(rnaturalearth)

source("scripts/00_plot_styles.R")

# Create output directories if they do not exist
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 0. Configuration & Parameters
# ==============================================================================
# Define global thresholds to ensure reproducibility and perfectly match 16S parameters
ABUNDANCE_THRESHOLD <- 0.0001  # 0.01% relative abundance for detection
MIN_SAMPLES_MAP     <- 10      # Minimum samples per country to display on map
MIN_SAMPLES_BAR     <- 10     # Minimum samples per region for bar chart
TARGET_TAXON        <- "o__Gastranaerophilales"
BAR_FILL_COLOR      <- "gray"

# ==============================================================================
# 1. Load Data and Extract Target Taxon Abundance
# ==============================================================================
cat("Loading Metalog processed data...\n")

# Use fast reading for the large species table
species_data <- fread("data/metalog_motus3_processed.tsv.gz", showProgress = FALSE) %>% 
  column_to_rownames(var = "rowname")

metadata <- fread("data/metadata_with_groups_processed.tsv.gz") %>% as.data.frame()
taxa_data <- fread("data/mOTUs_3.0.0_GTDB_tax.tsv.gz") %>% as.data.frame()

if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"

cat("Extracting Target Taxon...\n")
# Map column names to taxonomy to find Gastranaerophilales
mapping <- data.frame(original_col = colnames(species_data)) %>%
  mutate(id = str_extract(original_col, "(?<=\\[)[^\\[\\]]+(?=\\]$)")) %>%
  left_join(taxa_data, by = "id")

target_motus <- mapping %>% filter(order == TARGET_TAXON) %>% pull(original_col)

if(length(target_motus) == 0) stop("Gastranaerophilales not found in Metalog data!")

# Calculate total abundance for the target order per sample
target_abundance_df <- data.frame(
  sample_alias = rownames(species_data),
  Relative_Abundance = rowSums(species_data[, target_motus, drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
)

# Merge with metadata and apply threshold
df_analysis <- metadata %>%
  inner_join(target_abundance_df, by = "sample_alias") %>%
  mutate(Detected = Relative_Abundance > ABUNDANCE_THRESHOLD)

# ==============================================================================
# 2. Generate Map Data
# ==============================================================================
cat("Generating Map...\n")

map_stats <- df_analysis %>%
  group_by(country) %>%
  filter(n() >= MIN_SAMPLES_MAP) %>%
  summarise(
    TotalSamples = n(),
    PositiveSamples = sum(Detected, na.rm = TRUE),
    Prevalence = (PositiveSamples / TotalSamples) * 100,
    .groups = "drop"
  ) %>%
  # Fix country names to match rnaturalearth map data
  mutate(country_map = case_when(
    country == "United States" ~ "United States of America",
    country == "Democratic Republic of the Congo" ~ "Dem. Rep. Congo",
    country == "Republic of Congo" ~ "Congo",
    country == "Cote d'Ivoire" ~ "Côte d'Ivoire",
    country == "Central African Republic" ~ "Central African Rep.",
    TRUE ~ country
  ))

world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(name != "Antarctica")

p_map <- ggplot(data = world_map %>% left_join(map_stats, by = c("name" = "country_map"))) +
  geom_sf(aes(fill = Prevalence), color = "white", linewidth = 0.05) + 
  coord_sf(expand = FALSE, datum = NA) +
  scale_fill_gradient2(
    low = "#74ADD1",
    mid = "#FFFFCF",
    high = "#F46D43",
    midpoint = 50,
    na.value = "grey90",
    name = "Prevalence (%)", 
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100) 
  ) +
  labs(
    title = expression(paste(bold("Metagenome")))
  ) +
  theme_void(base_size = 6) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(hjust = 0.5, size = 6, face = "bold", margin = margin(b = 2)),
    legend.title = element_text(vjust = 0.8, face = "bold", size = 5),
    legend.text = element_text(size = 5),
    legend.key.width = unit(5, "mm"),
    legend.key.height = unit(1.5, "mm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    plot.margin = margin(2, 2, 2, 2)
  )

# ==============================================================================
# 3. Generate Bar Chart Data
# ==============================================================================
cat("Generating Bar Chart...\n")

bar_stats <- df_analysis %>%
  filter(!is.na(region), region != "") %>%
  group_by(region) %>%
  filter(n() >= MIN_SAMPLES_BAR) %>%
  summarise(
    TotalSamples = n(),
    PositiveSamples = sum(Detected, na.rm = TRUE),
    Prevalence = (PositiveSamples / TotalSamples) * 100,
    .groups = "drop"
  )

p_bar <- ggplot(bar_stats, aes(x = reorder(region, -Prevalence), y = Prevalence, fill = Prevalence)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#74ADD1",
    mid = "#FFFFCF",
    high = "#F46D43",
    midpoint = 50,
    limits = c(0, 100) 
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 100)) +
  labs(
    x = NULL, y = "Prevalence (%)"
  ) +
  theme_nature(base_size = 6) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 5),
    axis.text.y = element_text(color = "black", size = 5),
    axis.title.y = element_text(size = 5),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 6, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 10),
    legend.position = "none"
  )

# ==============================================================================
# 4. Export Plots and Tables
# ==============================================================================
cat("Exporting plots separately...\n")

outdir <- "results/figures/"

ggsave(paste0(outdir, "Figure3a_Metalog_Map.pdf"), p_map, 
       width = 60, height = 45, units = "mm", useDingbats = FALSE, bg = "transparent")
ggsave(paste0(outdir, "Figure3a_Metalog_Map.png"), p_map,
       width = 60, height = 45, units = "mm", dpi = 300, bg = "white")

ggsave(paste0(outdir, "Figure3a_Metalog_Region.pdf"), p_bar, 
       width = 35, height = 45, units = "mm", useDingbats = FALSE, bg = "transparent")
ggsave(paste0(outdir, "Figure3a_Metalog_Region.png"), p_bar,
       width = 35, height = 45, units = "mm", dpi = 300, bg = "white")

# Save Supplementary Tables
cat("Saving supplementary tables...\n")
write_tsv(df_analysis %>% select(sample_alias, country, region, Relative_Abundance, Detected), 
          "results/tables/Supplementary_Table_Metalog_Filtered_Samples.tsv")
write_tsv(map_stats, "results/tables/Supplementary_Table_Metalog_Country_Prevalence.tsv")
write_tsv(bar_stats, "results/tables/Supplementary_Table_Metalog_Regional_Prevalence.tsv")

cat("Analysis complete! Outputs saved to results/ directory.\n")

# ==============================================================================
# Generate Supplementary Tables
# ==============================================================================
cat("Saving supplementary tables...\n")

# --- Supplementary Table 6: Country-level Prevalence ---
supp_table_6 <- map_stats %>%
  select(
    Country = country,
    Total_Samples = TotalSamples,
    Positive_Samples = PositiveSamples,
    Prevalence_Percentage = Prevalence
  ) %>%
  mutate(Prevalence_Percentage = round(Prevalence_Percentage, 2)) %>%
  arrange(desc(Prevalence_Percentage), Country)

write_csv(supp_table_6, "results/tables/Supplementary_Table6_Country_Prevalence.csv")

cat("Analysis complete! Outputs saved to results/ directory.\n")
