# ==============================================================================
# Script: 07_lmm_higher_ranks.R
# Author: Kotaro Murai
# Description: 
#   Perform LMM Analysis on higher taxonomic ranks (Phylum to Genus).
#   Generate a consolidated Lollipop Plot (Figure 2c).
# ==============================================================================

# Load libraries
library(tidyverse)
library(data.table)
library(lmerTest)
library(broom.mixed)
library(future)
library(furrr)
library(emmeans)
library(ggnewscale) # For multiple color scales

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
    mutate(id = str_extract(original_col, "(?<=\\[).+?(?=\\])")) %>%
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
    Analysis_Group = factor(Analysis_Group, levels = c("Industrialized", "Non-Africa Non-Industrialized", "Africa")),
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
# Function to run LMER (Same logic as Script 05)
run_lmer_rank <- function(taxon_name, dataset) {
  subset_dt <- dataset[taxon == taxon_name]
  tryCatch({
    suppressMessages({
      model <- lmer(log10_abundance ~ Analysis_Group + disease_group + (1 | study) + (1 | country), data = subset_dt)
      emm <- emmeans(model, ~ Analysis_Group, lmer.df = "satterthwaite")
      contrast_results <- contrast(emm, method = list(AvsND = c(0, -1, 1), AvsNDev = c(-1, 0, 1)))
    })
    tidy(contrast_results) %>% mutate(taxon = taxon_name, .before = 1)
  }, error = function(e) return(NULL))
}

plan(multicore)
options(future.globals.maxSize = 8 * 1024^3)

# Storage for all results
all_ranks_results <- list()

for (rank in target_ranks) {
  cat(sprintf("\nProcessing Rank: %s ... \n", rank))
  
  # 1. Aggregate
  rank_df <- aggregate_by_rank(species_data, taxa_data, rank)
  
  # 2. Screening (Group-wise logic)
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
  
  # 3. Prepare Long Format
  long_dt <- bind_cols(model_metadata, rank_df %>% select(all_of(taxa_to_keep))) %>%
    pivot_longer(cols = all_of(taxa_to_keep), names_to = "taxon", values_to = "abundance") %>%
    mutate(log10_abundance = log10(abundance + 1E-6)) %>% 
    as.data.table()
  setkey(long_dt, taxon)
  
  # 4. Run LMM
  rank_results_list <- future_map(taxa_to_keep, ~run_lmer_rank(.x, long_dt), .progress = TRUE)
  
  # 5. Format & Store
  rank_res <- rbindlist(purrr::compact(rank_results_list)) %>%
    group_by(contrast) %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) %>%
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

plot_data <- df_all %>%
  filter(contrast == "AvsND") %>%
  filter(!str_detect(taxon, "Not_annotated|Incongruent|-1|Unassigned")) %>%
  left_join(master_tax_map, by = c("taxon" = "taxon_name")) %>%
  mutate(
    clean_name = str_remove(taxon, "^[pcofg]__"),
    Type = if_else(estimate > 0, "Enriched", "Depleted"),
    Rank = factor(Rank, levels = c("phylum", "class", "order", "family", "genus"))
  ) %>%
  filter(!str_detect(clean_phylum, "Unassigned|NA"))

top_n_show <- 15 
est_threshold <- 0.3

plot_enriched <- plot_data %>%
  filter(Type == "Enriched", q.value < 0.05, estimate > est_threshold) %>%
  arrange(desc(estimate)) %>% slice_head(n = top_n_show)

plot_depleted <- plot_data %>%
  filter(Type == "Depleted", q.value < 0.05, estimate < -est_threshold) %>%
  arrange(estimate) %>% slice_head(n = top_n_show)

final_plot_data <- bind_rows(plot_enriched, plot_depleted) %>%
  arrange(estimate) %>%
  distinct(taxon, .keep_all = TRUE) %>%
  mutate(
    unique_key = factor(taxon, levels = unique(taxon)),
    q_label = sprintf("FDR = %.1e", q.value)
  )

# Apply Master Colors
cols_status <- c("Enriched" = thesis_colors[["Africa"]], "Depleted" = thesis_colors[["Non-Africa Non-Industrialized"]])

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

rank_shapes <- c("phylum" = 22, "class" = 24, "order" = 23, "family" = 21, "genus" = 21)
rank_abbrev <- c("phylum"="(p)", "class"="(c)", "order"="(o)", "family"="(f)", "genus"="(g)")
y_labels <- setNames(
  paste0(final_plot_data$clean_name, " ", rank_abbrev[as.character(final_plot_data$Rank)]),
  final_plot_data$unique_key
)

p <- ggplot(final_plot_data, aes(x = estimate, y = unique_key)) +
  geom_vline(xintercept = 0, color = "gray50") +
  geom_vline(xintercept = c(est_threshold, -est_threshold), linetype = "dashed", color = "gray80") +
  
  geom_segment(aes(x = 0, xend = estimate, y = unique_key, yend = unique_key, color = Type), linewidth = 1) +
  scale_color_manual(values = cols_status, guide = "none") + 
  
  geom_text(aes(
    label = q_label, 
    x = estimate + sign(estimate) * 0.05, 
    hjust = ifelse(estimate > 0, 0, 1)    
  ), size = 3.5, color = "gray30") +
  
  geom_point(aes(fill = clean_phylum, shape = Rank), size = 4, color = "black", stroke = 0.5) + 
  scale_fill_manual(values = final_phylum_colors, name = "Phylum", guide = guide_legend(override.aes = list(shape = 21))) +
  scale_shape_manual(values = rank_shapes, name = "Rank") +
  
  scale_y_discrete(labels = y_labels) +
  scale_x_continuous(expand = expansion(mult = 0.25)) + 
  labs(
    title = "Top Differentially Abundant Taxa (All Ranks)",
    x = "LMM Estimate (Africa vs Non-Ind.)\n<-- Depleted in Africa       Enriched in Africa -->",
    y = NULL
  ) +
  
  theme_thesis(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "italic", size = 11, color = "black"),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray95")
  )

outfile <- "results/figures/Figure2c_Lollipop_AllRanks"
ggsave(paste0(outfile, ".pdf"), p, width = 12, height = 10, useDingbats = FALSE)
ggsave(paste0(outfile, ".png"), p, width = 12, height = 10, dpi = 300, bg = "white")