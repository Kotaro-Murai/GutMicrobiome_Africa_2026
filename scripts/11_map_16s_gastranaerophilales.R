
# ==============================================================================
# Script: 11_map_16s_gastranaerophilales.R
# Description: Plot world map of Gastranaerophilales prevalence using external 16S data.
#              (Figure 3d) - Streaming Version
# ==============================================================================

library(tidyverse)
library(data.table)
library(sf)
library(rnaturalearth)
library(countrycode)

# Load custom styles
source("scripts/00_plot_styles.R")
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load and Clean Metadata
# ==============================================================================
cat("Loading and cleaning metadata...\n")

pr <- fread("data/projects.csv.gz") %>% as.data.frame()
pr_clean <- pr %>%
  filter(
    !str_detect(subjects, regex("infant|mice|younger than 2|newborn|MICE|babies", ignore_case = TRUE)),
    avg_length >= 100,
    sample_type == "stool"
  )

md <- fread("data/sample_metadata.tsv.gz") %>% as.data.frame()
md_clean <- md %>%
  filter(
    project %in% pr_clean$project,
    iso != "unknown",
    iso != "",
    !is.na(iso)
  ) %>%
  mutate(id = paste0(project, "_", srr))

valid_ids <- unique(md_clean$id)
cat(sprintf("Metadata processed. Target samples: %d\n", length(valid_ids)))

# ==============================================================================
# 2. Streaming Process (Chunked Reading)
# ==============================================================================
cat("Starting streaming process (Reading file in chunks)...\n")
tax_file <- "data/taxonomic_table.csv.gz"

header <- names(fread(tax_file, nrows = 0))
target_taxon <- "Gastranaerophilales"
target_cols <- grep(target_taxon, header, value = TRUE)

if(length(target_cols) == 0) stop("Target taxon not found in header!")
cat(sprintf("Target taxon columns found: %d\n", length(target_cols)))

process_chunk <- function(chunk, pos) {
  # 🌟 ユーザー修正を反映: IDは2列目
  ids <- chunk[[2]]
  
  keep_idx <- which(ids %in% valid_ids)
  if (length(keep_idx) == 0) return(NULL) 
  
  chunk_sub <- chunk[keep_idx, ]
  ids_sub <- ids[keep_idx]
  
  if(length(target_cols) > 1) {
    target_val <- rowSums(chunk_sub[, target_cols, drop=FALSE], na.rm = TRUE)
  } else {
    target_val <- chunk_sub[[target_cols]]
  }
  
  # 🌟 ユーザー修正を反映: 総リード計算時は1,2列目を除外
  total_reads <- rowSums(chunk_sub[, -c(1,2), drop=FALSE], na.rm = TRUE)
  
  data.frame(
    id = ids_sub,
    Target_Count = target_val,
    Total_Reads = total_reads,
    stringsAsFactors = FALSE
  )
}

df_abundance <- suppressMessages(read_csv_chunked(
  tax_file,
  callback = DataFrameCallback$new(process_chunk),
  chunk_size = 10000,
  show_col_types = FALSE,
  progress = TRUE
))

if (is.null(df_abundance) || nrow(df_abundance) == 0) {
  stop("df_abundance is empty. Check if IDs match between metadata and taxonomic table.")
}

cat("Streaming complete. Calculating prevalence...\n")

# ==============================================================================
# 3. Normalization & Thresholding
# ==============================================================================
df_analysis <- md_clean %>%
  inner_join(df_abundance, by = "id") %>%
  mutate(
    Relative_Abundance = Target_Count / Total_Reads,
    Detected = Relative_Abundance > 0.0001 # 閾値 0.01%
  )

# ==============================================================================
# 4. Summarize by ISO Country Code
# ==============================================================================
cat("Aggregating by ISO code...\n")

# 16Sデータはサンプルが多いため、国ごとの信頼性を担保する
min_samples <- 50 

prevalence_stats <- df_analysis %>%
  group_by(iso) %>%
  filter(n() >= min_samples) %>%
  summarise(
    TotalSamples = n(),
    PositiveSamples = sum(Detected, na.rm = TRUE),
    Prevalence = (PositiveSamples / TotalSamples) * 100, # 地図用にパーセント化
    .groups = "drop"
  )

# ISOコードを地図データ標準の「3文字コード（iso3c）」に変換
prevalence_stats <- prevalence_stats %>%
  mutate(iso_a3 = countrycode(iso, origin = "iso2c", destination = "iso3c"))

# ==============================================================================
# 5. Visualization (World Map)
# ==============================================================================
cat("Plotting world map...\n")

world_map <- ne_countries(scale = "medium", returnclass = "sf")

map_data <- world_map %>% 
  left_join(prevalence_stats, by = "iso_a3")

p_map <- ggplot(data = map_data) +
  geom_sf(aes(fill = Prevalence), color = "white", linewidth = 0.2) +
  
  scale_fill_distiller(
    palette = "YlOrRd",
    direction = 1,
    na.value = "grey90",
    name = "Prevalence (%)",
    labels = function(x) paste0(x, "%")
  ) +
  
  labs(
    title = expression(paste(bold("Prevalence of "), bolditalic("Gastranaerophilales"), bold(" (16S Validation)"))),
    subtitle = paste0("External dataset | Threshold > 0.01% | Samples/Country >= ", min_samples)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40", margin = margin(b = 20)),
    legend.key.width = unit(2, "cm"),
    legend.title = element_text(vjust = 0.8, face = "bold")
  )

# Save
outfile <- "results/figures/Figure3d_Map_16S_Validation"
ggsave(paste0(outfile, ".pdf"), p_map, width = 12, height = 7, useDingbats = FALSE)
ggsave(paste0(outfile, ".png"), p_map, width = 12, height = 7, dpi = 300, bg = "white")

cat("Successfully saved Figure 3d (16S World Map)!\n")