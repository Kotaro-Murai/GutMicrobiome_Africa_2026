# ==============================================================================
# Script: Supple_phylogenetic_distribution.R
# Description: Phylogenetic distribution of specific genera across animal species.
#              Generates Supplementary Figures 5a, 5b, and 5c.
# ==============================================================================

# --- 0. Setup and Libraries ---
library(tidyverse)
library(data.table)
library(ggtree)
library(rotl)
library(ape)
library(patchwork)

# Load custom plotting styles
source("scripts/00_plot_styles.R")

# Create output directories
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# Global Parameters
MIN_SAMPLE_SIZE <- 10
ABUNDANCE_THRESHOLD <- 0.0001 # 0.01% relative abundance threshold for detection

# Targets for Supplementary Figure 3
target_list <- list(
  list(name = "Prevotella",    fig = "Supplementary_Figure5a_Prevotella"),
  list(name = "Succinivibrio", fig = "Supplementary_Figure5b_Succinivibrio"),
  list(name = "Treponema",     fig = "Supplementary_Figure5c_Treponema")
)

# ==============================================================================
# 1. Global Data Loading & Tree Construction (Run Once for Efficiency)
# ==============================================================================
cat("Loading animal and human datasets...\n")

# --- Animal Data ---
df_meta  <- read_tsv("data/metalog_animal_metadata.20251210.tsv.gz", show_col_types = FALSE)
df_motus <- fread("data/metalog_animal_motus3.20251210.tsv.gz") %>% as.data.frame()
df_tax   <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
colnames(df_tax)[1] <- "mOTU_id"

df_meta_clean <- df_meta %>%
  mutate(species_name = str_extract(tax_id_name, "^\\S+\\s+\\S+")) %>%
  filter(!is.na(species_name))

valid_species <- df_meta_clean %>%
  group_by(species_name) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= MIN_SAMPLE_SIZE) %>%
  pull(species_name)

# --- Fetch Phylogenetic Tree ---
cat("Fetching phylogenetic tree from Open Tree of Life...\n")
prevalence_base <- df_meta_clean %>% filter(species_name %in% valid_species) %>% distinct(species_name)
taxa_search <- tnrs_match_names(names = prevalence_base$species_name, context_name = "Animals")
valid_taxa <- taxa_search %>% filter(!is.na(ott_id))
tree <- tol_induced_subtree(ott_ids = valid_taxa$ott_id)

tree_clean <- tree
tree_clean$tip.label <- str_remove(tree_clean$tip.label, "_ott\\d+") %>% 
  str_replace_all("_", " ") %>%         
  str_remove(" \\(species in domain Eukaryota\\)")
tree_clean <- ladderize(tree_clean, right = FALSE)

# --- Human Data ---
human_species_data <- fread("data/metalog_motus3_processed.tsv.gz", showProgress = FALSE) %>% 
  column_to_rownames(var = "rowname")
human_metadata <- fread("data/metadata_with_groups_processed.tsv.gz") %>% as.data.frame()
human_taxa <- fread("data/mOTUs_3.0.0_GTDB_tax.tsv.gz") %>% as.data.frame()
if(!"id" %in% colnames(human_taxa)) colnames(human_taxa)[1] <- "id"

# ==============================================================================
# 2. Functionalize Plotting for Multiple Targets
# ==============================================================================

generate_phylogeny_plot <- function(target_genus, fig_name) {
  cat(sprintf("\nProcessing %s for %s...\n", target_genus, fig_name))
  
  # --- Animal Cohort Process ---
  target_ids_tax <- df_tax %>% filter(str_detect(genus, target_genus)) %>% pull(mOTU_id)
  clean_target_ids <- str_extract(as.character(target_ids_tax), "\\d+$")
  
  motus_cols <- colnames(df_motus)
  motus_numeric_ids <- str_extract(motus_cols, "\\d+(?=\\]$)")
  target_indices <- which(motus_numeric_ids %in% clean_target_ids)
  target_cols_names <- motus_cols[target_indices]
  
  total_counts <- rowSums(df_motus %>% select(-rowname), na.rm = TRUE)
  if(length(target_cols_names) > 0) {
    target_counts <- rowSums(df_motus %>% select(all_of(target_cols_names)), na.rm = TRUE)
  } else {
    target_counts <- rep(0, nrow(df_motus))
  }
  
  animal_abundance_df <- data.frame(
    sample_alias = df_motus$rowname,
    target_rel_abun = target_counts / total_counts
  ) %>% mutate(target_rel_abun = replace_na(target_rel_abun, 0)) 
  
  prevalence_data <- df_meta_clean %>%
    filter(species_name %in% valid_species) %>%
    inner_join(animal_abundance_df, by = "sample_alias") %>%
    group_by(species_name) %>%
    summarise(
      total_samples = n(),                    
      positive_samples = sum(target_rel_abun > ABUNDANCE_THRESHOLD),  
      mean_abundance = mean(target_rel_abun),
      .groups = "drop"
    ) %>%
    mutate(
      prevalence = positive_samples / total_samples, 
      tax_id_name = species_name 
    )
  
  # --- Human Cohort Process ---
  human_mapping <- data.frame(original_col = colnames(human_species_data)) %>%
    mutate(id = str_extract(original_col, "(?<=\\[)[^\\[\\]]+(?=\\]$)")) %>%
    left_join(human_taxa, by = "id")
  
  human_target_motus <- human_mapping %>% 
    filter(str_detect(genus, target_genus)) %>% 
    pull(original_col)
  
  if(length(human_target_motus) > 0) {
    human_abundance_df <- data.frame(
      sample_alias = rownames(human_species_data),
      Relative_Abundance = rowSums(human_species_data[, human_target_motus, drop = FALSE], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  } else {
    human_abundance_df <- data.frame(sample_alias = rownames(human_species_data), Relative_Abundance = 0)
  }
  
  human_df <- human_metadata %>%
    inner_join(human_abundance_df, by = "sample_alias") %>%
    mutate(Detected = Relative_Abundance > ABUNDANCE_THRESHOLD) %>%
    filter(!is.na(Analysis_Group), Analysis_Group %in% c("Sub-Saharan Africa", "Other Non-Industrialized", "Industrialized")) %>%
    group_by(Analysis_Group) %>%
    summarise(
      total_samples = n(),
      positive_samples = sum(Detected, na.rm = TRUE),
      prevalence = positive_samples / total_samples,
      .groups = "drop"
    ) %>%
    mutate(
      Analysis_Group = factor(Analysis_Group, levels = c("Sub-Saharan Africa", "Other Non-Industrialized", "Industrialized")),
      lbl = paste0("italic('Homo sapiens')~'(", Analysis_Group, ", ", total_samples, ")'")
    ) %>%
    arrange(Analysis_Group)
  
  # --- Visualization ---
  prev_limits <- c(0, 1)
  prev_breaks <- seq(0, 1, by = 0.25)
  size_range  <- c(0.8, 2.5) 
  
  p <- ggtree(tree_clean, layout = "rectangular", linewidth = 0.25) %<+% prevalence_data
  tree_max_x <- max(p$data$x, na.rm = TRUE)
  plot_max_x <- tree_max_x * 2
  
  p_tree <- p +
    geom_tree(linewidth = 0.25, color = "gray40") +
    geom_tiplab(
      aes(label = paste0("italic('", tax_id_name, "')~'(", total_samples, ")'")),
      parse = TRUE, align = TRUE, linesize = 0.2, linetype = "dotted", 
      offset = tree_max_x * 0.05, size = 2
    ) +
    geom_tippoint(
      aes(size = prevalence, fill = prevalence), 
      shape = 21, color = "gray20", stroke = 0.25
    ) +
    scale_fill_gradient2(
      low = "#74ADD1", mid = "#FFFFCF", high = "#F46D43", midpoint = 0.5, 
      na.value = "grey90", name = "Prevalence", breaks = prev_breaks,
      labels = scales::percent_format(accuracy = 1), limits = c(0, 1),
      guide = guide_legend(reverse = TRUE) 
    ) +
    scale_size_continuous(
      range = size_range, name = "Prevalence", breaks = prev_breaks, 
      labels = scales::percent_format(accuracy = 1), limits = c(0, 1),
      guide = guide_legend(reverse = TRUE) 
    ) +
    xlim(0, plot_max_x) + 
    theme_nature(base_size = 6) + 
    theme(
      axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      panel.border = element_blank(), panel.grid.major = element_blank(),
      legend.position = "bottom", legend.box = "horizontal",
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6), legend.key.size = unit(3, "mm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
      plot.margin = margin(t = 2, r = 2, b = 0, l = 2, unit = "mm")
    )
  
  p_human <- ggplot(human_df, aes(x = 0, y = fct_rev(Analysis_Group))) + 
    geom_segment(aes(x = 0, xend = tree_max_x * 0.5, y = Analysis_Group, yend = Analysis_Group), 
                 linetype = "dotted", color = "gray40", linewidth = 0.2) +
    geom_point(
      aes(size = prevalence, fill = prevalence), 
      shape = 21, color = "gray20", stroke = 0.25
    ) +
    geom_text(
      aes(x = tree_max_x * 0.6, label = lbl), parse = TRUE, 
      hjust = 0, size = 2 
    ) +
    scale_fill_gradient2(
      low = "#74ADD1", mid = "#FFFFCF", high = "#F46D43", midpoint = 0.5,
      na.value = "grey90", limits = c(0, 1)
    ) +
    scale_size_continuous(range = size_range, limits = prev_limits) +
    xlim(0, plot_max_x) +
    theme_nature(base_size = 6) + 
    theme(
      axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), panel.border = element_blank(),
      panel.grid.major = element_blank(), legend.position = "none",
      plot.margin = margin(t = 0, r = 2, b = 2, l = 2, unit = "mm") 
    )
  
  final_plot <- (p_tree / p_human) + plot_layout(heights = c(15, 2))
  
  # --- Export ---
  outfile <- paste0("results/figures/", fig_name)
  ggsave(paste0(outfile, ".pdf"), final_plot, width = 90, height = 120, units = "mm", useDingbats = FALSE)
  ggsave(paste0(outfile, ".png"), final_plot, width = 90, height = 120, units = "mm", dpi = 300, bg = "white")
  
  # Export Summary Table for this Genus
  sample_abundance <- animal_abundance_df %>%
    mutate(Detected = target_rel_abun > ABUNDANCE_THRESHOLD) %>%
    rename(!!paste0(target_genus, "_Abundance") := target_rel_abun)
  
  supp_table <- df_meta_clean %>%
    filter(species_name %in% valid_species) %>%
    left_join(sample_abundance, by = "sample_alias") %>%
    select(Sample_ID = sample_alias, Study = study, Species_Name = species_name, Age = age, Diet = diet, 
           ends_with("_Abundance"), Detected) %>%
    mutate(across(ends_with("_Abundance"), ~ signif(., 4))) %>%
    arrange(Species_Name, Study, desc(Detected), Sample_ID)
  
  write_csv(supp_table, paste0("results/tables/Supplementary_Table_", target_genus, ".csv"))
  cat(sprintf("Successfully saved %s!\n", fig_name))
}

# ==============================================================================
# 3. Execute Loop over Targets
# ==============================================================================

for (target in target_list) {
  generate_phylogeny_plot(target_genus = target$name, fig_name = target$fig)
}

# Clean up memory
rm(human_species_data, human_metadata, human_taxa)
gc()

cat("\nAll Supplementary Figures generated successfully!\n")