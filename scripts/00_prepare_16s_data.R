# ==============================================================================
# Script: 00_prepare_16s_data.R
# Description: Filter raw metadata and taxonomic count table to define the 16S 
#              cohort. Converts raw counts to relative abundance and exports 
#              compressed processed data for downstream analysis.
# ==============================================================================

library(tidyverse)
library(data.table)

cat("Starting 16S cohort metadata preprocessing...\n")

# ==============================================================================
# 1. Metadata Filtering
# ==============================================================================
pr <- fread("data/projects.csv.gz") %>% as.data.frame()
md <- fread("data/sample_metadata.tsv.gz") %>% as.data.frame()

# Exclusion criteria: remove infants, animal models, and short reads
pr_clean <- pr %>%
  filter(
    !str_detect(subjects, regex("infant|mice|younger than 2|newborn|MICE|babies", ignore_case = TRUE)),
    avg_length >= 100,
    sample_type == "stool"
  )

md_clean <- md %>%
  filter(
    project %in% pr_clean$project,
    iso != "unknown",
    iso != "",
    !is.na(iso)
  ) %>%
  mutate(iso = if_else(iso == "UM", "US", iso),
         id = paste0(project, "_", srr))

valid_ids <- unique(md_clean$id)

write_tsv(md_clean, "data/16s_metadata_processed.tsv.gz")

cat(sprintf("Metadata processing complete. Retained %d valid samples.\n", length(valid_ids)))

# ==============================================================================
# 2. Taxonomic Table Filtering & Relative Abundance Conversion (Streaming)
# ==============================================================================
cat("Filtering taxonomic table and calculating relative abundance...\n")

tax_file_in <- "data/taxonomic_table.csv.gz"
tax_file_out <- "data/16s_taxonomic_table_processed.csv.gz"

if (file.exists(tax_file_out)) file.remove(tax_file_out)

process_and_write_chunk <- function(chunk, pos) {
  ids <- chunk[[2]] 
  keep_idx <- which(ids %in% valid_ids)
  
  if (length(keep_idx) > 0) {
    chunk_sub <- chunk[keep_idx, ]
    
    meta_cols <- chunk_sub[, 1:2]
    count_cols <- chunk_sub[, -c(1, 2)]
    
    total_reads <- rowSums(count_cols, na.rm = TRUE)
    total_reads[total_reads == 0] <- 1 
    
    rel_abund <- sweep(count_cols, 1, total_reads, "/")
    
    chunk_out <- bind_cols(meta_cols, rel_abund)
    
    write_csv(
      chunk_out, 
      file = tax_file_out, 
      append = TRUE, 
      col_names = (pos == 1)
    )
  }
}

suppressMessages(read_csv_chunked(
  tax_file_in,
  callback = SideEffectChunkCallback$new(process_and_write_chunk),
  chunk_size = 10000,
  show_col_types = FALSE,
  progress = TRUE
))

cat("Preparation complete! Processed files saved in 'data/'\n")