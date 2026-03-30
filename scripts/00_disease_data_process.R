# ==============================================================================
# Script: 00_disease_data_process.R
# Author: Kotaro Murai
# Description: 
#   Clean and standardize the 'disease' column in the metadata.
#   Groups synonymous or highly related diseases into broader categories
#   and handles NA values as 'Control'.
# ==============================================================================

# Load libraries
library(tidyverse)

# --- 1. Load Data ---
cat("Loading raw metadata...\n")

# Set file paths
input_file  <- "data/metalog_metadata_processed.tsv.gz"
output_file <- "data/metalog_metadata_processed.tsv.gz" 

# Load metadata
md <- read_tsv(input_file, show_col_types = FALSE)

# --- 2. Process Disease Data ---
cat("Standardizing disease labels...\n")

# Define disease groupings
md_cleaned <- md %>%
  mutate(
    
    # Standardize and group diseases
    disease = case_when(
      disease == "Cardiovascular system disease" ~ "Cardiovascular disease",
      
      disease %in% c("Metabolic dysfunction-associated steatohepatitis", 
                     "Metabolic associated fatty liver disease", 
                     "Vascular liver", 
                     "Metabolic dysfunction-associated steatotic liver disease", 
                     "Liver disease", 
                     "Metabolic liver disease") ~ "Chronic liver disease",
      
      disease %in% c("Adenocarcinoma of the lung", 
                     "Non-small cell lung cancer") ~ "Lung cancer",
      
      disease == "Impaired glucose tolerance" ~ "Prediabetes",
      
      disease == "End-stage renal disease" ~ "Chronic kidney disease",
      
      # Keep other diseases as they are
      TRUE ~ disease
    )
  )

# --- 3. Save Output ---
cat("Saving cleaned metadata...\n")

# Write back to compressed TSV
write_tsv(md_cleaned, output_file)

# Print a quick summary to verify the changes
cat("\n--- Disease count summary (Top 15) ---\n")
md_cleaned %>% 
  count(disease, sort = TRUE) %>% 
  slice_head(n = 15) %>% 
  print(n = Inf)

cat("\nSuccessfully finished Script 00_disease_data_process.R!\n")