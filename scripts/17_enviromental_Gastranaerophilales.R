# ==============================================================================
# Script: 17_environmental_absence.R
# Description: Demonstrates the global absence of the target order 
#              (Gastranaerophilales) across diverse environmental and marine 
#              datasets. Figure 5b.
# Author: Kotaro Murai
# ==============================================================================

# --- 0. Setup and Libraries ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(sf)
  library(rnaturalearth)
  library(patchwork)
})

# Load custom plotting styles
source("scripts/00_plot_styles.R")

# Create output directories
outdir <- "results/figures/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Global Parameters
TARGET_ORDER <- "Gastranaerophilales"
ABUNDANCE_THRESHOLD <- 0.0001 # 0.01% relative abundance threshold for detection

# --- Taxonomy File ---
cat("Loading taxonomy database...\n")
tax_file <- "data/mOTUs_3.0.0_GTDB_tax.tsv.gz"
df_tax <- read_tsv(tax_file, show_col_types = FALSE)
if(!"id" %in% colnames(df_tax)) colnames(df_tax)[1] <- "id"


# ==============================================================================
# 1. Functions for Data Extraction
# ==============================================================================

# Function to safely process large mOTUs tables and calculate target abundance
extract_target_abundance <- function(meta_path, motus_path, dataset_name) {
  cat(sprintf(" -> Processing %s dataset...\n", dataset_name))
  
  # Load metadata
  meta <- read_tsv(meta_path, show_col_types = FALSE)
  
  # Ensure 'environment_material' exists to prevent mapping errors
  if(!"environment_material" %in% colnames(meta)) {
    meta <- meta %>% mutate(environment_material = NA_character_)
  }
  
  # Extract essential columns and enforce data types
  meta <- meta %>%
    select(sample_alias, longitude, latitude, environment_material) %>%
    mutate(
      sample_alias = as.character(sample_alias),
      longitude = as.numeric(longitude),
      latitude = as.numeric(latitude),
      environment_material = as.character(environment_material)
    )
  
  # Load mOTUs profile
  motus <- fread(motus_path) %>% as.data.frame() %>% column_to_rownames("rowname")
  
  # Map columns to taxonomy and identify target columns
  mapping <- data.frame(original_col = colnames(motus)) %>%
    mutate(id = str_extract(original_col, "(?<=\\[).+?(?=\\])")) %>%
    left_join(df_tax, by = "id")
  
  target_cols <- mapping %>% filter(str_detect(order, TARGET_ORDER)) %>% pull(original_col)
  
  total_reads <- rowSums(motus, na.rm = TRUE)
  target_reads <- rowSums(motus[, target_cols, drop = FALSE], na.rm = TRUE)
  
  abundance_df <- data.frame(
    sample_alias = rownames(motus),
    Abundance_Target = target_reads / total_reads,
    Dataset = dataset_name,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Abundance_Target = replace_na(Abundance_Target, 0),
      Detected_Target = Abundance_Target > ABUNDANCE_THRESHOLD 
    )
  
  # Merge with metadata
  res <- meta %>% inner_join(abundance_df, by = "sample_alias")
  
  # Free memory
  rm(motus)
  gc()
  
  return(res)
}

# ==============================================================================
# 2. Data Integration and Biome Categorization
# ==============================================================================
cat("\nIntegrating environmental and oceanic datasets...\n")

# File paths
env_meta_file  <- "data/metalog_environmental_metadata.20251210.tsv.gz"
env_motus_file <- "data/metalog_environmental_motus3.20251210.tsv.gz"
ocean_meta_file  <- "data/metalog_ocean_metadata.20251210.tsv.gz"
ocean_motus_file <- "data/metalog_ocean_motus3.20251210.tsv.gz"

# Extract data
df_env   <- extract_target_abundance(env_meta_file, env_motus_file, "Environment")
df_ocean <- extract_target_abundance(ocean_meta_file, ocean_motus_file, "Ocean")

# Combine datasets and apply rigorous biome categorization based on ENVO ontology
df_analysis <- bind_rows(df_env, df_ocean) %>%
  mutate(
    # Clean string for regex matching
    mat_clean = tolower(str_remove(environment_material, "\\s*\\[ENVO:.*\\]"))
  ) %>%
  mutate(
    # Categorization logic representing major global biomes
    Biome_Group = case_when(
      # Ocean dataset partitioning
      Dataset == "Ocean" & str_detect(mat_clean, "sediment") ~ "Sediment", 
      Dataset == "Ocean" ~ "Ocean",                                      
      
      # Environmental dataset mapping
      str_detect(mat_clean, "soil") ~ "Soil",                                # [ENVO:00001998]
      str_detect(mat_clean, "sediment") ~ "Sediment",                        # [ENVO:00002007]
      str_detect(mat_clean, "fresh water|river|lake|groundwater") ~ "Freshwater", # [ENVO:00002011]
      str_detect(mat_clean, "waste|sewage|sludge") ~ "Wastewater",           # [ENVO:00002001]
      str_detect(mat_clean, "air|dust|aerosol") ~ "Air & Dust",              # [ENVO:00002005], [ENVO:00002008]
      
      TRUE ~ "Other"
    )
  )

cat(sprintf("Total global samples processed: %d\n", nrow(df_analysis)))

# ==============================================================================
# 3. Visualization Panel A: Global Sampling Maps (6 Biomes)
# ==============================================================================
cat("\nGenerating Global Maps for Specific Biomes (Panel A)...\n")

# Fetch global map and remove Antarctica
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin != "Antarctica") 

# Filter and correct mapping coordinates
df_map <- df_analysis %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  filter(longitude >= -180 & longitude <= 180) %>%
  filter(latitude >= -90 & latitude <= 90) %>%
  arrange(Detected_Target) 

target_biomes <- c("Ocean", "Soil", "Freshwater", "Wastewater", "Air & Dust", "Sediment")

generate_and_save_map <- function(biome_name) {
  df_sub <- df_map %>% filter(Biome_Group == biome_name)
  
  safe_name <- str_replace_all(biome_name, "[^A-Za-z0-9]", "_") %>% str_replace_all("__+", "_") %>% str_remove_all("^_|_$")
  
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "gray90", color = "white", linewidth = 0.1) +
    geom_point(data = filter(df_sub, !Detected_Target), 
               aes(x = longitude, y = latitude), 
               color = "lightblue", alpha = 0.9, size = 0.2, stroke = 0) +
    geom_point(data = filter(df_sub, Detected_Target), 
               aes(x = longitude, y = latitude), 
               color = "#D55E00", alpha = 0.9, size = 0.3, stroke = 0) +
    coord_sf(expand = FALSE, xlim = c(-180, 180), ylim = c(-80, 90)) + 
    labs(title = biome_name) +
    theme_void() +
    theme(
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5, margin = margin(b = 1, unit = "mm")),
      plot.margin = margin(1, 1, 1, 1, unit = "mm")
    )
  
  outfile <- paste0(outdir, "Figure5b_Map_", safe_name)
  ggsave(paste0(outfile, ".pdf"), p_map, width = 45, height = 30, units = "mm", useDingbats = FALSE)
  ggsave(paste0(outfile, ".png"), p_map, width = 45, height = 30, units = "mm", dpi = 300, bg = "white")
  
  cat(sprintf(" -> Saved %s Map (45x30mm).\n", biome_name))
}
walk(target_biomes, generate_and_save_map)


# ==============================================================================
# 4. Visualization Panel B: Prevalence Bar Chart (Horizontal)
# ==============================================================================
cat("Generating Biome Prevalence Chart (Panel B)...\n")

# Calculate prevalence metrics
bar_stats <- df_analysis %>%
  group_by(Biome_Group) %>%
  summarise(
    n_samples = n(),
    Prevalence = sum(Detected_Target) / n(),
    .groups = "drop"
  ) %>%
  mutate(
    Biome_Label = paste0(Biome_Group, "(n=", n_samples, ")")
  )

# Extract levels for ordering
levels_ordered <- bar_stats %>%
  filter(Biome_Group != "Other") %>%
  arrange(desc(n_samples)) %>%
  pull(Biome_Label)

level_other <- bar_stats %>%
  filter(Biome_Group == "Other") %>%
  pull(Biome_Label)

bar_stats <- bar_stats %>%
  mutate(Biome_Label = factor(Biome_Label, levels = rev(c(levels_ordered, level_other))))

max_prev <- max(bar_stats$Prevalence) * 100
y_limit <- if(max_prev < 1) 1 else ceiling(max_prev)

biome_colors <- c(
  "Ocean"      = "#1f78b4", 
  "Soil"       = "#33a02c", 
  "Freshwater" = "#a6cee3",
  "Wastewater" = "#fb9a99", 
  "Air & Dust" = "#fdbf6f", 
  "Sediment"   = "#b2df8a", 
  "Other"      = "gray70"  
)
p_bar <- ggplot(bar_stats, aes(x = Prevalence * 100, y = Biome_Label, fill = Biome_Group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = biome_colors, guide = "none") + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, y_limit)) +
  labs(x = "Prevalence (%)", y = NULL) +
  theme_nature(base_size = 6) +
  theme(
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title.x = element_text(margin = margin(t = 2, unit = "mm")),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm")
  )

# ==============================================================================
# 5. Export Bar Plot
# ==============================================================================
cat("\nExporting Bar Chart...\n")

bar_outfile <- paste0(outdir, "Figure5b_Bar_Environmental_Gastraranerophilales")

ggsave(paste0(bar_outfile, ".pdf"), p_bar, width = 90, height = 30, units = "mm", useDingbats = FALSE)
ggsave(paste0(bar_outfile, ".png"), p_bar, width = 90, height = 30, units = "mm", dpi = 300, bg = "white")
cat(" -> Saved Bar Chart (90x30mm).\n")


# ==============================================================================
# 6. Export Supplementary Table 13 (Environmental Samples Metadata & Detection)
# ==============================================================================
cat("\nExporting Supplementary Table 13 (Environmental Samples Metadata)...\n")
dir.create("results/tables/", showWarnings = FALSE, recursive = TRUE)

supp_table_13 <- df_analysis %>%
  select(
    Sample_ID = sample_alias,
    Dataset,
    Assigned_Biome_Group = Biome_Group,
    Original_Environment_Material = environment_material,
    Longitude = longitude,
    Latitude = latitude,
    Gastranaerophilales_Abundance = Abundance_Target,
    Detected = Detected_Target
  ) %>%
  mutate(
    Longitude = round(Longitude, 4),
    Latitude = round(Latitude, 4),
    Gastranaerophilales_Abundance = signif(Gastranaerophilales_Abundance, 4)
  ) %>%
  arrange(Assigned_Biome_Group, Dataset, desc(Detected), Sample_ID)

write_csv(supp_table_13, "results/tables/Supplementary_Table_13_Environmental_Samples.csv")
cat(" -> Saved Supplementary Table 13 (Sample-level Metadata).\n")

cat("\nAnalysis pipeline complete. All outputs saved to results/ directories.\n")