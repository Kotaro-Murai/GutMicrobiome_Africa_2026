# ==============================================================================
# Script: 05_lmm_species_volcano.R
# Author: Kotaro Murai
# Description: 
#   Perform Linear Mixed Models (LMM) on species-level relative abundances.
#   Identify Africa-specific enriched/depleted species and generate a Volcano Plot.
# ==============================================================================

# Load libraries
library(tidyverse)
library(data.table)
library(lmerTest)
library(broom.mixed)
library(future)
library(furrr)
library(emmeans)
library(ggrepel)

# Load custom plot styles
source("scripts/00_plot_styles.R")

# Create output directories
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PART 1: Data Preparation & LMM Analysis
# ==============================================================================

cat("Loading data...\n")
Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

# Load abundance table
species_data <- read_tsv("data/metalog_motus3_processed.tsv.gz", show_col_types = FALSE) %>%
  column_to_rownames(var = "rowname")

# Load pre-grouped metadata from Step 02 (Group defined as Analysis_Group)
metadata <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)

# Load taxonomy
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"
if("species" %in% colnames(taxa_data)) taxa_data <- taxa_data %>% rename(taxon_species = species)

# --- 1. Metadata Formatting for LMM ---
cat("Preparing metadata for LMER...\n")

model_metadata <- metadata %>%
  select(sample_alias, country, study, disease, Analysis_Group) %>%
  mutate(
    # Set Factor levels for Analysis_Group (Important for contrasts)
    # Order: 1=Industrialized, 2=Non-Ind, 3=Africa
    Analysis_Group = factor(Analysis_Group, levels = c("Industrialized", "Non-Africa Non-Industrialized", "Africa")),
    
    # Handle NA diseases as Control
    disease = replace_na(disease, "Control")
  )

# Group rare diseases into 'Other' to keep all samples
disease_counts <- model_metadata %>% count(disease)
major_diseases <- disease_counts %>% filter(n >= 100 | disease == "Control") %>% pull(disease)

model_metadata <- model_metadata %>%
  mutate(
    disease_group = if_else(disease %in% major_diseases, disease, "Other"),
    disease_group = factor(disease_group),
    disease_group = relevel(disease_group, ref = "Control") # Control is the baseline
  )

cat(sprintf("Success! Kept all %d samples.\n", nrow(model_metadata)))

# --- 2. Species Screening ---
cat("Screening major species per group (Prevalence > 10%, Mean Abundance > 0.01%)...\n")

prevalence_threshold <- 0.1
abundance_threshold <- 0.0001
species_to_keep <- c()

for (g in levels(model_metadata$Analysis_Group)) {
  idx <- which(model_metadata$Analysis_Group == g)
  group_data <- species_data[idx, , drop = FALSE]
  prev <- colMeans(group_data > 0, na.rm = TRUE)
  abun <- colMeans(group_data, na.rm = TRUE)
  keep_in_group <- names(which(prev >= prevalence_threshold & abun >= abundance_threshold))
  species_to_keep <- unique(c(species_to_keep, keep_in_group))
  
  cat(sprintf("  - %s: %d species met the threshold.\n", g, length(keep_in_group)))
}

cat("Species screening finished. Total remaining species across all groups:", length(species_to_keep), "\n")

# --- 3. LMM Execution Function ---
run_lmer_and_emmeans <- function(species_name, dataset) {
  subset_dt <- dataset[species == species_name]
  
  result <- tryCatch({
    suppressMessages({
      # 1. LMM Model (Accounting for study and country as random effects)
      model <- lmer(log10_abundance ~ Analysis_Group + disease_group + (1 | study) + (1 | country), data = subset_dt)
      
      # 2. Estimated Marginal Means
      emm <- emmeans(model, ~ Analysis_Group, lmer.df = "satterthwaite")
      
      # 3. Contrasts
      # Based on Factor levels: 1:Ind, 2:Non-Ind, 3:Africa
      contrast_results <- contrast(emm, 
                                   method = list(
                                     AvsND   = c(0, -1, 1), # Africa vs Non-Industrialized
                                     AvsNDev = c(-1, 0, 1)  # Africa vs Industrialized
                                   ))
    })
    
    # 4. Tidy results
    tidy(contrast_results) %>% mutate(species = species_name, .before = 1)
    
  }, error = function(e) { return(NULL) })
  return(result)
}

# --- 4. Main LMM Loop (Parallelized) ---
cat("Starting final analysis loop with emmeans...\n")

# Prepare data in long format (data.table for speed)
long_dt <- bind_cols(model_metadata, species_data %>% select(all_of(species_to_keep))) %>%
  pivot_longer(cols = all_of(species_to_keep), names_to = "species", values_to = "abundance") %>%
  mutate(log10_abundance = log10(abundance + 1E-6)) %>%
  as.data.table()

setkey(long_dt, species)

plan(multicore) 
options(future.globals.maxSize = 40 * 1024^3)

# Run LMM for all kept species
lmm_results_list <- future_map(species_to_keep, ~run_lmer_and_emmeans(.x, long_dt), .progress = TRUE)

plan(sequential)
rm(long_dt, species_data); gc()

# --- 5. Format and Save Results ---
cat("Aggregating LMM results...\n")
all_contrasts <- rbindlist(purrr::compact(lmm_results_list))

# Save full results
write_csv(all_contrasts, "results/tables/lmm_results_species_full.csv")

# Apply FDR correction (BH method) and pivot wider
results_wide <- all_contrasts %>%
  select(species, contrast, estimate, p.value) %>%
  group_by(contrast) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  pivot_wider(id_cols = species, 
              names_from = contrast, 
              values_from = c(estimate, p.value, q.value))

write_csv(results_wide, "results/tables/lmm_results_species_wide.csv")

# Extract Africa-specific species (Compared to BOTH groups)
africa_specific_enriched <- results_wide %>%
  filter(q.value_AvsND < 0.05 & estimate_AvsND > 0.3 &
           q.value_AvsNDev < 0.05 & estimate_AvsNDev > 0.3) %>%
  arrange(desc(estimate_AvsND))

africa_specific_depleted <- results_wide %>%
  filter(q.value_AvsND < 0.05 & estimate_AvsND < -0.3 &
           q.value_AvsNDev < 0.05 & estimate_AvsNDev < -0.3) %>%
  arrange(estimate_AvsND)

# Save lists for downstream phylogenetic tree (Script 06)
write_csv(africa_specific_enriched, "results/tables/africa_enriched_species.csv")
write_csv(africa_specific_depleted, "results/tables/africa_depleted_species.csv")


# ==============================================================================
# PART 2: Volcano Plot (Figure 2a)
# ==============================================================================

cat("Generating Volcano Plot...\n")

# Prepare plotting data (Focusing on Africa vs Non-Industrialized 'AvsND')
plot_data <- results_wide %>%
  mutate(mOTU_id = str_extract(species, "(?<=\\[).+?(?=\\])")) %>%
  left_join(taxa_data, by = c("mOTU_id" = "id")) %>%
  select(
    species,
    estimate = estimate_AvsND, 
    q_value = q.value_AvsND,
    taxon_species
  ) %>%
  mutate(
    log10_q = -log10(q_value),
    clean_name = coalesce(taxon_species, species) %>% 
      str_remove("\\[.*\\]") %>% 
      str_remove("^s__") %>%
      str_trim()
  )

# Thresholds
q_value_threshold <- 0.05
estimate_threshold <- 0.3 

plot_data <- plot_data %>%
  mutate(
    significance = case_when(
      q_value < q_value_threshold & estimate > estimate_threshold  ~ "Enriched in Africa",
      q_value < q_value_threshold & estimate < -estimate_threshold ~ "Depleted in Africa",
      TRUE                                                         ~ "Not significant"
    ),
    # Only label highly significant features to avoid clutter
    label_text = if_else(
      q_value < 0.05 & abs(estimate) > 0.3, 
      clean_name,                              
      NA_character_                            
    )
  )

# Plotting
volcano_plot <- ggplot(plot_data, aes(x = estimate, y = log10_q, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  
  geom_text_repel(aes(label = label_text),
                  max.overlaps = 15,
                  size = 3,
                  color = "black",
                  box.padding = 0.5,
                  point.padding = 0.2) +
  
  geom_vline(xintercept = c(-estimate_threshold, estimate_threshold), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(q_value_threshold), linetype = "dashed", color = "gray50") +
  
  scale_color_manual(values = c(
    "Enriched in Africa" = thesis_colors[["Africa"]],  # Match Africa color from 00_plot_styles
    "Depleted in Africa" = thesis_colors[["Non-Africa Non-Industrialized"]],  # Match Non-Ind color
    "Not significant"    = "grey80"
  ), name = "Significance") +
  
  labs(
    title = "Volcano Plot: Africa vs. Non-Africa Non-Industrialized",
    subtitle = "Linear Mixed Model Analysis (Random effect: Country & Study)",
    x = "LMM Estimate (log10 scale)\n<-- Africa Depleted      Africa Enriched -->",
    y = expression("-Log"[10]*"(FDR)")
  ) +
  
  theme_thesis(base_size = 14) +
  theme(
    legend.position = "top",
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title.x = element_text(size = 12)
  )

# Display and save
print(volcano_plot)
ggsave("results/figures/Figure2a_Volcano_AvsND.pdf", plot = volcano_plot, width = 8, height = 7, useDingbats = FALSE)
ggsave("results/figures/Figure2a_Volcano_AvsND.png", plot = volcano_plot, width = 8, height = 7, dpi = 300)

cat("Successfully finished Script 05!\n")