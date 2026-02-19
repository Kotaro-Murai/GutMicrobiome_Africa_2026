# ==============================================================================
# Script: 09_validation_16s_prevalence.R
# Description: Validation of Gastranaerophilales prevalence using external 16S data.
#              (Figure 3b) - Streaming Processing Version
# ==============================================================================

library(tidyverse)
library(data.table)

# Load custom styles
source("scripts/00_plot_styles.R")
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load and Clean Metadata
# ==============================================================================
cat("Loading and cleaning metadata...\n")

# Load Projects
pr <- fread("data/projects.csv.gz") %>% as.data.frame()

# Filter Projects
pr_clean <- pr %>%
  filter(
    !str_detect(subjects, regex("infant|mice|younger than 2|newborn|MICE|babies", ignore_case = TRUE)),
    avg_length >= 100,
    sample_type == "stool"
  )

# Load Sample Metadata
md <- fread("data/sample_metadata.tsv.gz") %>% as.data.frame()

# Filter Samples
md_clean <- md %>%
  filter(
    project %in% pr_clean$project,
    iso != "unknown",
    !is.na(region),
    region != ""
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

cat("Streaming complete. Calculating prevalence...\n")

# ==============================================================================
# 3. Normalization & Thresholding
# ==============================================================================

df_analysis <- md_clean %>%
  inner_join(df_abundance, by = "id")

df_analysis <- df_analysis %>%
  mutate(
    Relative_Abundance = Target_Count / Total_Reads,
    Detected = Relative_Abundance > 0.0001 
  )

# ==============================================================================
# 4. Summarize by Region
# ==============================================================================
prevalence_stats <- df_analysis %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  filter(n() >= 100) %>%
  summarise(
    TotalSamples = n(),
    PositiveSamples = sum(Detected),
    Prevalence = PositiveSamples / TotalSamples,
    .groups = "drop"
  ) %>%
  arrange(desc(Prevalence))

print(head(prevalence_stats))

# ==============================================================================
# 5. Visualization (Neutral Style)
# ==============================================================================
cat("Plotting...\n")

neutral_color <- "steelblue"

p <- ggplot(prevalence_stats, aes(x = reorder(region, -Prevalence), y = Prevalence)) +
  
  geom_bar(stat = "identity", width = 0.8, fill = neutral_color) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.1))) +
  
  labs(
    title = expression(paste("Detection Frequency of ", bolditalic("Gastranaerophilales"), " (16S Validation)")),
    subtitle = paste0("External dataset (N=", sum(prevalence_stats$TotalSamples), "), Threshold > 0.01%"),
    x = NULL,
    y = "Detection Rate (%)"
  ) +
  
  theme_thesis(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

# Save
outfile <- "results/figures/Figure3b_Validation_16S_Prevalence"
ggsave(paste0(outfile, ".pdf"), p, width = 8, height = 6, useDingbats = FALSE)
ggsave(paste0(outfile, ".png"), p, width = 8, height = 6, dpi = 300, bg = "white")

cat("Successfully saved Figure 3b\n")