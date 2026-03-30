# ==============================================================================
# Script: 10_16s_Map_Gastranaerophilales.R
# Description: Generates map and bar chart of Gastranaerophilales prevalence 
#              (Figure 3a) using pre-processed 16S rRNA sequencing data.
# Author: Kotaro Murai
# ==============================================================================

library(tidyverse)
library(data.table)
library(sf)
library(rnaturalearth)
library(countrycode)
library(openxlsx) 

# Load custom plotting themes
source("scripts/00_plot_styles.R")

# Create output directories if they do not exist
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 0. Configuration & Parameters
# ==============================================================================
# Define global thresholds to ensure reproducibility and easy modification
ABUNDANCE_THRESHOLD <- 0.0001  # 0.01% relative abundance for detection
MIN_SAMPLES_MAP     <- 10      # Minimum samples per country to display on map
MIN_SAMPLES_BAR     <- 10      # Minimum samples per region for bar chart
BAR_FILL_COLOR      <- "gray"

# ==============================================================================
# 1. Load Data and Extract Target Taxon Abundance
# ==============================================================================
cat("Loading pre-processed data...\n")

md_clean <- fread("data/16s_metadata_processed.tsv.gz") %>% as.data.frame()
tax_file <- "data/16s_taxonomic_table_processed.csv.gz"

header <- names(fread(tax_file, nrows = 0))
target_taxon <- "Gastranaerophilales"
target_cols <- grep(target_taxon, header, value = TRUE)

if(length(target_cols) == 0) stop("Target taxon not found in header!")

# Function to process taxonomic table in chunks to preserve memory
process_chunk <- function(chunk, pos) {
  if(length(target_cols) > 1) {
    target_val <- rowSums(chunk[, target_cols, drop = FALSE], na.rm = TRUE)
  } else {
    target_val <- chunk[[target_cols]]
  }
  data.frame(
    id = chunk[[2]],
    Relative_Abundance = target_val,
    stringsAsFactors = FALSE
  )
}

df_abundance <- suppressMessages(read_csv_chunked(
  tax_file, callback = DataFrameCallback$new(process_chunk),
  chunk_size = 10000, show_col_types = FALSE
))

df_analysis <- md_clean %>%
  inner_join(df_abundance, by = "id") %>%
  mutate(
    Detected = Relative_Abundance > ABUNDANCE_THRESHOLD,
    region_metalog = countrycode(iso, origin = "iso2c", destination = "region") 
  )

# ==============================================================================
# 2. Generate Map Data (Panel A)
# ==============================================================================
cat("Generating Map...\n")

map_stats <- df_analysis %>%
  group_by(iso) %>%
  filter(n() >= MIN_SAMPLES_MAP) %>%
  summarise(
    TotalSamples = n(),
    PositiveSamples = sum(Detected, na.rm = TRUE),
    Prevalence = (PositiveSamples / TotalSamples) * 100,
    .groups = "drop"
  ) %>%
  mutate(iso_a3 = countrycode(iso, origin = "iso2c", destination = "iso3c"))

# Map data with patch for France and Norway ISO codes
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(name != "Antarctica") %>%
  mutate(iso_a3 = case_when(
    admin == "France" ~ "FRA",
    admin == "Norway" ~ "NOR",
    TRUE ~ iso_a3
  ))

p_map <- ggplot(data = world_map %>% left_join(map_stats, by = "iso_a3")) +
  geom_sf(aes(fill = Prevalence), color = "white", linewidth = 0.1) + 
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
    title = expression(paste(bold("16S rRNA")))
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
# 3. Generate Bar Chart Data (Panel B)
# ==============================================================================
cat("Generating Bar Chart...\n")

bar_stats <- df_analysis %>%
  filter(!is.na(region_metalog), region_metalog != "") %>%
  group_by(region_metalog) %>%
  filter(n() >= MIN_SAMPLES_BAR) %>%
  summarise(
    TotalSamples = n(),
    PositiveSamples = sum(Detected, na.rm = TRUE),
    Prevalence = PositiveSamples / TotalSamples,
    .groups = "drop"
  )

p_bar <- ggplot(bar_stats, aes(x = reorder(region_metalog, -Prevalence), y = Prevalence)) +
  geom_bar(stat = "identity", width = 0.7, fill = BAR_FILL_COLOR) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL, y = "Prevalence (%)"
  ) +
  theme_nature(base_size = 6) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 5),
    axis.text.y = element_text(color = "black", size = 5),
    axis.title.y = element_text(size = 5),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(2, 2, 2, 2)
  )

# ==============================================================================
# 4. Export Plots and Tables
# ==============================================================================
cat("Exporting plots separately...\n")

outdir <- "results/figures/"

ggsave(paste0(outdir, "Figure3a_16S_Map.pdf"), p_map, 
       width = 70, height = 45, units = "mm", useDingbats = FALSE, bg = "transparent")
ggsave(paste0(outdir, "Figure3a_16S_Map.png"), p_map,
       width = 70, height = 45, units = "mm", dpi = 300, bg = "white")

ggsave(paste0(outdir, "Figure3a_16S_Region.pdf"), p_bar, 
       width = 20, height = 45, units = "mm", useDingbats = FALSE, bg = "transparent")
ggsave(paste0(outdir, "Figure3a_16S_Region.png"), p_bar,
       width = 20, height = 45, units = "mm", dpi = 300, bg = "white")

# Save Supplementary Tables
cat("Saving supplementary tables...\n")
write_tsv(df_analysis %>% select(id, iso, region_metalog, Relative_Abundance, Detected), 
          "results/tables/Supplementary_Table_16S_Filtered_Samples.tsv")
write_tsv(map_stats, "results/tables/Supplementary_Table_16S_Country_Prevalence.tsv")
write_tsv(bar_stats, "results/tables/Supplementary_Table_16S_Regional_Prevalence.tsv")

cat("Analysis complete! Outputs saved to results/ directory.\n")


# ==============================================================================
# Generate Supplementary Table 4 (Combined Metagenome & 16S)
# ==============================================================================
cat("Saving combined supplementary table...\n")

supp_table_16s <- map_stats %>%
  mutate(Country = countrycode(iso, origin = "iso2c", destination = "country.name")) %>%
  select(
    Country,
    Total_Samples = TotalSamples,
    Positive_Samples = PositiveSamples,
    Prevalence_Percentage = Prevalence
  ) %>%
  mutate(Prevalence_Percentage = round(Prevalence_Percentage, 2)) %>%
  arrange(desc(Prevalence_Percentage), Country)

metagenome_file <- "results/tables/Supplementary_Table_4_Country_Prevalence.csv"
if(file.exists(metagenome_file)) {
  supp_table_metagenome <- read_csv(metagenome_file, show_col_types = FALSE)
} else {
  warning("Metagenome data not found. Only 16S data will be saved.")
  supp_table_metagenome <- data.frame(Message = "Run Script 09 first to generate this sheet.")
}

sheets <- list(
  "Metagenome" = supp_table_metagenome,
  "16S_rRNA" = supp_table_16s
)

write.xlsx(sheets, "results/tables/Supplementary_Table_4_Gastranaerophilales_Prevalence.xlsx")
