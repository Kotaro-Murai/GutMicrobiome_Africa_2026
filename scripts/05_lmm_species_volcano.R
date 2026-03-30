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
library(patchwork)

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
    # Order: 1=Industrialized, 2=Non-Ind, 3=Sub-Saharan Africa
    Analysis_Group = factor(Analysis_Group, levels = c("Industrialized", "Other Non-Industrialized", "Sub-Saharan Africa")),
    
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
      # 1. LMM Model
      model <- lmer(log10_abundance ~ Analysis_Group + disease_group + (1 | study) + (1 | country), data = subset_dt)
      
      # 2. Estimated Marginal Means
      emm <- emmeans(model, ~ Analysis_Group, lmer.df = "satterthwaite")
      
      # 3. Contrasts
      contrast_results <- contrast(emm, 
                                   method = list(
                                     SSAvsONI = c(0, -1, 1),
                                     SSAvsInd = c(-1, 0, 1)
                                   ))
      res_contrast <- tidy(contrast_results) %>% 
        select(contrast, estimate, p.value)
      res_emm <- tidy(emm) %>% 
        rename(contrast = Analysis_Group) %>% 
        mutate(contrast = paste0("EMM_", contrast)) %>% 
        select(contrast, estimate) 
      bind_rows(res_contrast, res_emm) %>% 
        mutate(species = species_name, .before = 1)
    })
  }, error = function(e) { return(NULL) })
  return(result)
}

# --- 4. Main LMM Loop (Parallelized) ---
cat("Starting final analysis loop with emmeans...\n")

# Prepare data in long format (data.table for speed)
long_dt <- bind_cols(model_metadata, species_data %>% select(all_of(species_to_keep))) %>%
  pivot_longer(cols = all_of(species_to_keep), names_to = "species", values_to = "abundance") %>%
  mutate(log10_abundance = log10(abundance + 1E-4)) %>%
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
all_contrasts <- rbindlist(purrr::compact(lmm_results_list), fill = TRUE)

# Save full results
write_csv(all_contrasts, "results/tables/lmm_results_species_full.csv")

# Apply FDR correction (BH method) and pivot wider
results_wide <- all_contrasts %>%
  select(species, contrast, estimate, p.value) %>%
  group_by(contrast) %>%
  mutate(q.value = if_else(!is.na(p.value), p.adjust(p.value, method = "BH"), NA_real_)) %>%
  ungroup() %>%
  pivot_wider(id_cols = species, 
              names_from = contrast, 
              values_from = c(estimate, p.value, q.value))

write_csv(results_wide, "results/tables/lmm_results_species_wide.csv")

# ==============================================================================
# 6. Prepare Supplementary Table 2 (Full LMM Results with Taxonomy)
# ==============================================================================
cat("Generating Supplementary Table 2...\n")

supp_table_2 <- results_wide %>%
  mutate(mOTU_id = str_extract(species, "(?<=\\[).+?(?=\\])")) %>%
  left_join(taxa_data, by = c("mOTU_id" = "id")) %>%
  select(
    mOTU_id, species, 
    domain, phylum, class, order, family, genus, taxon_species,
    # SSA vs ONI
    estimate_SSAvsONI, p.value_SSAvsONI, q.value_SSAvsONI,
    # SSA vs Ind
    estimate_SSAvsInd, p.value_SSAvsInd, q.value_SSAvsInd,
    # EMM (Estimated Marginal Means)
    contains("EMM")
  ) %>%
  rename(
    Estimate_SSA_vs_ONI = estimate_SSAvsONI,
    Pvalue_SSA_vs_ONI = p.value_SSAvsONI,
    FDR_SSA_vs_ONI = q.value_SSAvsONI,
    Estimate_SSA_vs_Ind = estimate_SSAvsInd,
    Pvalue_SSA_vs_Ind = p.value_SSAvsInd,
    FDR_SSA_vs_Ind = q.value_SSAvsInd
  ) %>%
  arrange(FDR_SSA_vs_ONI)

write_csv(supp_table_2, "results/tables/Supplementary_Table_2_LMM_results.csv")
cat("Supplementary Table 2 saved.\n")


# ==============================================================================
# PART 2: Volcano Plots (Figure 2a) - Stacked (SSA vs ONI & SSA vs Ind)
# ==============================================================================

cat("Generating Stacked Volcano Plots...\n")
library(patchwork)

results_wide <- read_csv("results/tables/lmm_results_species_wide.csv", show_col_types = FALSE)

# 1. Thresholds
q_value_threshold <- 0.05
estimate_threshold <- 0.1 

# 2. Prepare plotting data (Base data without strict dual-filtering for colors)
plot_data <- results_wide %>%
  mutate(mOTU_id = str_extract(species, "(?<=\\[).+?(?=\\])")) %>%
  left_join(taxa_data, by = c("mOTU_id" = "id"))

# 3. Colors
volcano_colors <- c(
  "Enriched" = LMM_colors[["Enriched"]],
  "Depleted" = LMM_colors[["Depleted"]],
  "Not significant" = "gray85" 
)

# 4. Create Base Plot Function
create_volcano <- function(df, x_col, y_col, x_label) {
  
  df_plot <- df %>%
    mutate(
      significance = case_when(
        .data[[y_col]] < q_value_threshold & .data[[x_col]] > estimate_threshold ~ "Enriched",
        .data[[y_col]] < q_value_threshold & .data[[x_col]] < -estimate_threshold ~ "Depleted",
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(significance != "Not significant")
  
  ggplot(df_plot, aes(x = .data[[x_col]], y = -log10(.data[[y_col]]), color = significance)) +
    geom_vline(xintercept = c(-estimate_threshold, estimate_threshold), linetype = "dashed", color = "gray50", linewidth = 0.25) +
    geom_hline(yintercept = -log10(q_value_threshold), linetype = "dashed", color = "gray50", linewidth = 0.25) +
    geom_point(alpha = 0.8, size = 0.8, stroke = 0) +
    
    scale_color_manual(values = volcano_colors) +
    
    labs(x = x_label, y = expression("-Log"[10]*"(FDR)")) +
    
    theme_nature(base_size = 6) +
    theme(
      legend.position = "none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# 5. Generate Individual Plots
p_oni <- create_volcano(
  plot_data, 
  x_col = "estimate_SSAvsONI", 
  y_col = "q.value_SSAvsONI", 
  x_label = "LMM estimate\n(Sub-Saharan Africa vs. Other Non-Industrialized)"
)

p_ind <- create_volcano(
  plot_data, 
  x_col = "estimate_SSAvsInd", 
  y_col = "q.value_SSAvsInd", 
  x_label = "LMM estimate\n(Sub-Saharan Africa vs. Industrialized)"
)

# 6. Combine with Patchwork
combined_volcano <- (p_oni / p_ind) + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.text = element_text(size = 7),
    legend.margin = margin(t = 0, b = 0)
  )

# 7. Save
ggsave("results/figures/Figure2a_Volcano.pdf", plot = combined_volcano, width = 70, height = 125, units = "mm", useDingbats = FALSE)
ggsave("results/figures/Figure2a_Volcano.png", plot = combined_volcano, width = 70, height = 125, units = "mm", dpi = 300)

cat("Successfully finished Script 05!\n")