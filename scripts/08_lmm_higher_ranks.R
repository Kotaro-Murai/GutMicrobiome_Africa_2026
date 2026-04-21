# ==============================================================================
# Script: 08_lmm_higher_ranks.R
# Author: Kotaro Murai
# Description: 
#   Perform LMM Analysis on higher taxonomic ranks (Phylum to Genus).
# ==============================================================================

# Load libraries
library(tidyverse)
library(data.table)
library(lmerTest)
library(broom.mixed)
library(future)
library(furrr)
library(emmeans)

# Load custom plot styles
source("scripts/00_plot_styles.R")

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PART 1: Data Loading & Aggregation
# ==============================================================================
cat("Loading data...\n")
Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

# Load metadata with groups
metadata <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)

# Load species abundance (mOTUs3)
species_data <- read_tsv("data/metalog_motus3_processed.tsv.gz", show_col_types = FALSE) %>%
  column_to_rownames(var = "rowname")

# Load taxonomy
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"
if("species" %in% colnames(taxa_data)) taxa_data <- taxa_data %>% rename(taxon_species = species)

# Define ranks to analyze
target_ranks <- c("phylum", "class", "order", "family", "genus")

# --- Function to Aggregate Abundance by Rank ---
aggregate_by_rank <- function(df, tax_df, rank_level) {
  mapping <- data.frame(original_col = colnames(df)) %>%
    mutate(id = str_extract(original_col, "(?<=\\[)[^\\[\\]]+(?=\\][^\\[]*$)")) %>%
    left_join(tax_df, by = "id") %>%
    filter(!is.na(!!sym(rank_level)))
  
  rank_list <- split(mapping$original_col, mapping[[rank_level]])
  
  agg_list <- lapply(names(rank_list), function(r) {
    cols <- rank_list[[r]]
    if (length(cols) == 1) {
      res <- df[[cols]]
    } else {
      res <- rowSums(df[, cols], na.rm = TRUE)
    }
    setNames(data.frame(res), r)
  })
  
  agg_df <- bind_cols(agg_list)
  rownames(agg_df) <- rownames(df)
  return(agg_df)
}

# --- Prepare Metadata for Model ---
model_metadata <- metadata %>%
  select(sample_alias, country, study, disease, Analysis_Group) %>%
  mutate(
    Analysis_Group = factor(Analysis_Group, levels = c("Industrialized", "Other Non-Industrialized", "Sub-Saharan Africa")),
    disease = replace_na(disease, "Control")
  )

disease_counts <- model_metadata %>% count(disease)
major_diseases <- disease_counts %>% filter(n >= 100 | disease == "Control") %>% pull(disease)

model_metadata <- model_metadata %>%
  mutate(
    disease_group = if_else(disease %in% major_diseases, disease, "Other"),
    disease_group = factor(disease_group),
    disease_group = relevel(disease_group, ref = "Control")
  )

# ==============================================================================
# PART 2: LMM Analysis Loop
# ==============================================================================
run_lmer_rank <- function(taxon_name, dataset) {
  subset_dt <- dataset[taxon == taxon_name]
  tryCatch({
    suppressMessages({
      model <- lmer(log10_abundance ~ Analysis_Group + disease_group + (1 | study) + (1 | country), data = subset_dt)
      emm <- emmeans(model, ~ Analysis_Group, lmer.df = "satterthwaite")
      
      contrast_results <- contrast(emm, method = list(
        SSAvsONI = c(0, -1, 1), 
        SSAvsInd = c(-1, 0, 1)
      ))
      
      res_contrast <- tidy(contrast_results) %>% select(contrast, estimate, p.value)
      res_emm <- tidy(emm) %>% 
        rename(contrast = Analysis_Group) %>% 
        mutate(contrast = paste0("EMM_", contrast)) %>% 
        select(contrast, estimate)
      
      bind_rows(res_contrast, res_emm) %>% mutate(taxon = taxon_name, .before = 1)
      
    })
  }, error = function(e) return(NULL))
}

plan(multicore)
options(future.globals.maxSize = 8 * 1024^3)

all_ranks_results <- list()

for (rank in target_ranks) {
  cat(sprintf("\nProcessing Rank: %s ... \n", rank))
  
  rank_df <- aggregate_by_rank(species_data, taxa_data, rank)
  
  prevalence_threshold <- 0.1
  abundance_threshold <- 0.0001
  taxa_to_keep <- c()
  
  for (g in levels(model_metadata$Analysis_Group)) {
    idx <- which(model_metadata$Analysis_Group == g)
    group_data <- rank_df[idx, , drop = FALSE]
    prev <- colMeans(group_data > 0, na.rm = TRUE)
    abun <- colMeans(group_data, na.rm = TRUE)
    keep_in_group <- names(which(prev >= prevalence_threshold & abun >= abundance_threshold))
    taxa_to_keep <- unique(c(taxa_to_keep, keep_in_group))
  }
  cat(sprintf("  -> Retained %d taxa after screening.\n", length(taxa_to_keep)))
  
  if(length(taxa_to_keep) == 0) next
  
  long_dt <- bind_cols(model_metadata, rank_df %>% select(all_of(taxa_to_keep))) %>%
    pivot_longer(cols = all_of(taxa_to_keep), names_to = "taxon", values_to = "abundance") %>%
    mutate(log10_abundance = log10(abundance + 1E-4)) %>% 
    as.data.table()
  setkey(long_dt, taxon)
  
  rank_results_list <- future_map(taxa_to_keep, ~run_lmer_rank(.x, long_dt), .progress = TRUE)
  rank_res <- rbindlist(purrr::compact(rank_results_list), fill = TRUE) %>%
    group_by(contrast) %>%
    mutate(q.value = if_else(!is.na(p.value), p.adjust(p.value, method = "BH"), NA_real_)) %>%
    ungroup() %>%
    mutate(Rank = rank)
  
  write_csv(rank_res, sprintf("results/tables/lmm_results_%s.csv", rank))
  all_ranks_results[[rank]] <- rank_res
  
  rm(rank_df, long_dt); gc()
}

plan(sequential)

# ==============================================================================
# PART 3: Visualization (Consolidated Lollipop Plot)
# ==============================================================================
cat("\nGenerating Figure 2d (All Ranks Lollipop)...\n")

target_ranks <- c("phylum", "class", "order", "family", "genus")
df_all <- map_dfr(target_ranks, ~read_csv(sprintf("results/tables/lmm_results_%s.csv", .x), show_col_types = FALSE))

taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"

master_tax_map <- taxa_data %>%
  select(phylum, class, order, family, genus) %>%
  pivot_longer(cols = c(class, order, family, genus), names_to = "rank_type", values_to = "taxon_name") %>%
  select(taxon_name, phylum) %>%
  distinct() %>%
  bind_rows(taxa_data %>% distinct(phylum) %>% mutate(taxon_name = phylum)) %>%
  filter(!is.na(taxon_name)) %>%
  distinct(taxon_name, .keep_all = TRUE) %>%
  mutate(clean_phylum = str_remove(phylum, "^p__"))

# ==============================================================================
# Generate Supplementary Table 4 (Full LMM Results for Higher Ranks)
# ==============================================================================
cat("Generating Supplementary Table 4 (Higher Ranks LMM)...\n")

supp_table_4 <- df_all %>%
  pivot_wider(
    id_cols = c(Rank, taxon),
    names_from = contrast,
    values_from = c(estimate, p.value, q.value)
  ) %>%
  select(
    Rank, Taxon = taxon,
    # SSA vs ONI
    Estimate_SSA_vs_ONI = estimate_SSAvsONI,
    Pvalue_SSA_vs_ONI = p.value_SSAvsONI,
    FDR_SSA_vs_ONI = q.value_SSAvsONI,
    # SSA vs Ind
    Estimate_SSA_vs_Ind = estimate_SSAvsInd,
    Pvalue_SSA_vs_Ind = p.value_SSAvsInd,
    FDR_SSA_vs_Ind = q.value_SSAvsInd,
    # EMM (Estimated Marginal Means)
    contains("EMM")
  ) %>%
  mutate(Rank = factor(Rank, levels = c("phylum", "class", "order", "family", "genus"))) %>%
  arrange(Rank, FDR_SSA_vs_ONI)

write_csv(supp_table_4, "results/tables/Supplementary_Table4_HigherRanks_LMM.csv")
cat("Supplementary Table 4 saved.\n\n")

plot_data <- df_all %>%
  filter(contrast %in% c("SSAvsInd", "SSAvsONI")) %>%
  filter(!str_detect(taxon, "Not_annotated|Incongruent|-1|Unassigned")) %>%
  pivot_wider(
    id_cols = c(taxon, Rank),
    names_from = contrast,
    values_from = c(estimate, q.value)
  ) %>%
  filter(sign(estimate_SSAvsInd) == sign(estimate_SSAvsONI)) %>%
  mutate(
    conservative_est = if_else(
      abs(estimate_SSAvsInd) < abs(estimate_SSAvsONI),
      estimate_SSAvsInd,
      estimate_SSAvsONI
    ),
    Type = if_else(conservative_est > 0, "Enriched", "Depleted"),
    estimate = conservative_est
  ) %>%
  left_join(master_tax_map, by = c("taxon" = "taxon_name")) %>%
  mutate(
    clean_name = str_remove(taxon, "^[pcofg]__"),
    Rank = factor(Rank, levels = c("phylum", "class", "order", "family", "genus"))
  ) %>%
  filter(!str_detect(clean_phylum, "Unassigned|NA"))

top_n_show <- 10 
est_threshold <- 0.1

plot_enriched <- plot_data %>%
  filter(Type == "Enriched", q.value_SSAvsInd < 0.05 & q.value_SSAvsONI < 0.05, estimate > est_threshold) %>%
  arrange(desc(estimate)) %>% slice_head(n = top_n_show)

plot_depleted <- plot_data %>%
  filter(Type == "Depleted", q.value_SSAvsInd < 0.05 & q.value_SSAvsONI < 0.05, estimate < -est_threshold) %>%
  arrange(estimate) %>% slice_head(n = top_n_show)

final_plot_data <- bind_rows(plot_enriched, plot_depleted) %>%
  arrange(estimate) %>%
  distinct(taxon, .keep_all = TRUE) %>%
  mutate(unique_key = factor(taxon, levels = unique(taxon)))

# Apply Master Colors
cols_status <- c("Enriched" = LMM_colors[["Enriched"]], 
                 "Depleted" = LMM_colors[["Depleted"]])

used_phyla <- unique(final_plot_data$clean_phylum)
defined_cols <- phylum_master_colors[names(phylum_master_colors) %in% used_phyla]
missing_phyla <- setdiff(used_phyla, names(phylum_master_colors))
if(length(missing_phyla) > 0) {
  extra_cols <- taxa_colors_20[1:length(missing_phyla)]
  names(extra_cols) <- missing_phyla
  final_phylum_colors <- c(defined_cols, extra_cols)
} else {
  final_phylum_colors <- defined_cols
}

rank_shapes <- c("phylum" = 25, "class" = 24, "order" = 23, "family" = 22, "genus" = 21)
rank_abbrev <- c("phylum"="(p)", "class"="(c)", "order"="(o)", "family"="(f)", "genus"="(g)")
y_labels <- setNames(
  paste0(final_plot_data$clean_name, " ", rank_abbrev[as.character(final_plot_data$Rank)]),
  final_plot_data$unique_key
)

p <- ggplot(final_plot_data, aes(x = estimate, y = unique_key)) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.25) +
  geom_vline(xintercept = c(est_threshold, -est_threshold), linetype = "dashed", color = "gray80", linewidth = 0.25) +
  
  geom_segment(aes(x = 0, xend = estimate, y = unique_key, yend = unique_key, color = Type), linewidth = 0.5) +
  scale_color_manual(values = cols_status, guide = "none") + 
  
  geom_point(aes(fill = clean_phylum, shape = Rank), size = 2, color = "black", stroke = 0.25) + 
  scale_fill_manual(values = final_phylum_colors, guide = "none") +
  scale_shape_manual(values = rank_shapes, name = "Rank") +
  
  scale_y_discrete(labels = y_labels) +
  scale_x_continuous(expand = expansion(mult = 0.15)) + 
  labs(x = "LMM Estimate", y = NULL, title = NULL) +
  
  theme_nature(base_size = 6) +
  theme(
    axis.text.y = element_text(face = "italic", size = 7, color = "black"),
    
    # Embedded Legend styling
    legend.position = c(0.98, 0.02),        
    legend.justification = c(1, 0),         
    legend.background = element_blank(),    
    legend.key = element_blank(),           
    legend.box.background = element_blank(),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(4, "mm"),
    
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.25)
  )

outfile <- "results/figures/Figure2d_Lollipop_AllRanks"
ggsave(paste0(outfile, ".pdf"), p, width = 88, height = 80, units = "mm", useDingbats = FALSE, bg = "transparent")
ggsave(paste0(outfile, ".png"), p, width = 88, height = 80, units = "mm", dpi = 300, bg = "white")

cat("Successfully saved Figure 2d!\n")