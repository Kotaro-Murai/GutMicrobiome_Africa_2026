# ==============================================================================
# Script: 07_phylogenetic_tree.R
# Author: Kotaro Murai
# Description: 
#   Map LMM results onto the GTDB phylogenetic tree.
#   Features: Phylum ring and LMM Status ring. 
# ==============================================================================

# Load libraries
library(tidyverse)
library(ape)
library(ggtree)
library(treeio)
library(ggtreeExtra)
library(ggnewscale)
library(RColorBrewer)

# Load custom plot styles
source("scripts/00_plot_styles.R")

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load Data
# ==============================================================================
cat("Loading analysis results, taxonomy, and tree...\n")

# LMM Results from Script 05
res_df <- read_csv("results/tables/lmm_results_species_wide.csv", show_col_types = FALSE)

# Taxonomy data
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"
if("species" %in% colnames(taxa_data)) taxa_data <- taxa_data %>% rename(taxon_species = species)

# GTDB Tree & Metadata
tree_file <- "data/bac120_r207.tree"
meta_file <- "data/bac120_metadata_r207.tsv.gz"
if(!file.exists(tree_file)) stop("GTDB tree file not found in 'data/'!")

big_tree <- read.tree(tree_file)
gtdb_meta <- read_tsv(meta_file, col_select = c("accession", "gtdb_taxonomy", "gtdb_representative"), show_col_types = FALSE)

# ==============================================================================
# 2. Data Matching & Tree Pruning
# ==============================================================================
cat("Matching mOTUs to GTDB tree...\n")

# Clean GTDB metadata to get 's__Species_name' mapping
gtdb_map <- gtdb_meta %>%
  mutate(species_gtdb = str_extract(gtdb_taxonomy, "s__[^;]+$")) %>%
  arrange(species_gtdb, desc(gtdb_representative)) %>%
  distinct(species_gtdb, .keep_all = TRUE) %>%
  select(accession, species_gtdb)

# 1. Define LMM Status (Enriched / Depleted / NS)
res_clean <- res_df %>%
  filter(species != "-1") %>% 
  mutate(
    motu_id = str_extract(species, "(?<=\\[).+?(?=\\]$)"),
    Status = case_when(
      q.value_SSAvsInd < 0.05 & estimate_SSAvsInd > 0.1 & q.value_SSAvsONI < 0.05 & estimate_SSAvsONI > 0.1 ~ "Enriched",
      q.value_SSAvsInd < 0.05 & estimate_SSAvsInd < -0.1 & q.value_SSAvsONI < 0.05 & estimate_SSAvsONI < -0.1 ~ "Depleted",
      TRUE ~ "NS"
    )
  )

# 2. Match with taxa_data to get GTDB names
matched_step1 <- res_clean %>%
  inner_join(taxa_data, by = c("motu_id" = "id")) %>%
  mutate(
    species_key = if_else(str_detect(taxon_species, "^s__"), taxon_species, paste0("s__", taxon_species))
  )

# 3. Match with GTDB Tree Accessions
matched_data <- matched_step1 %>%
  inner_join(gtdb_map, by = c("species_key" = "species_gtdb")) %>%
  filter(accession %in% big_tree$tip.label)

if(nrow(matched_data) == 0) stop("No species matched the tree. Check GTDB version or ID formats.")

# 4. Prune Tree
pruned_tree <- keep.tip(big_tree, matched_data$accession)
cat("Successfully pruned tree. Remaining tips:", length(pruned_tree$tip.label), "\n")

# ==============================================================================
# 3. Prepare Visualization Data
# ==============================================================================
cat("Preparing tree metadata...\n")

viz_data <- matched_data %>%
  select(accession, species, species_key, phylum, Status) %>%
  mutate(
    phylum_group = str_remove(phylum, "^p__"),
    phylum_group = replace_na(phylum_group, "Unassigned"),
    Status = factor(Status, levels = c("Enriched", "Depleted", "NS")),
    clean_label = if_else(Status != "NS", 
                          str_replace_all(str_remove(species_key, "^s__"), "_", " "), 
                          NA_character_)
  )

d_tree <- data.frame(label = pruned_tree$tip.label) %>%
  left_join(viz_data, by = c("label" = "accession"))

# ==============================================================================
# 4. Plotting (Phylum & Status Rings)
# ==============================================================================
cat("Generating circular tree...\n")

# Only 2 rings now
global_x_levels <- c("Phylum", "Status")

d_tree <- d_tree %>%
  mutate(
    phylum_x = factor("Phylum", levels = global_x_levels),
    status_x = factor("Status", levels = global_x_levels)
  )

p_base <- ggtree(pruned_tree, layout = "fan", linewidth = 0.05, open.angle = 15) %>%
  rotate_tree(90)

# --- Color Definitions ---
phylum_counts <- table(d_tree$phylum_group)
sorted_phyla <- sort(unique(d_tree$phylum_group))
d_tree$phylum_group <- factor(d_tree$phylum_group, levels = sorted_phyla)

used_phyla <- unique(d_tree$phylum_group)
defined_cols <- phylum_master_colors[names(phylum_master_colors) %in% used_phyla]
missing_phyla <- setdiff(used_phyla, names(phylum_master_colors))
if(length(missing_phyla) > 0) {
  extra_cols <- taxa_colors_20[1:length(missing_phyla)]
  names(extra_cols) <- missing_phyla
  cols_phylum <- c(defined_cols, extra_cols)
} else {
  cols_phylum <- defined_cols
}
cols_phylum <- cols_phylum[sorted_phyla]

cols_status <- c(
  "Enriched" = LMM_colors[["Enriched"]], 
  "Depleted" = LMM_colors[["Depleted"]],
  "NS"       = "gray90"
)

# --- Build the Plot Layers ---
p <- p_base %<+% d_tree

# [Ring 1] Phylum
p1 <- p +
  geom_fruit(geom = geom_tile, 
             mapping = aes(fill = phylum_group, x = phylum_x), 
             pwidth = 0.2, 
             offset = 0.02,
             color = NA, 
             axis.params = list(
               axis = "x",
               text.size = 1.8,
               hjust = 0,
               vjust = 0.5,
               line.color = NA)
  ) +
  scale_fill_manual(
    values = cols_phylum, 
    name = "Phylum", 
    na.value = "gray95",
    guide = guide_legend(direction = "horizontal", 
                         nrow = 3, 
                         title.position = "top", 
                         order = 1) 
  ) +
  new_scale_fill()

# [Ring 2] Significance Status
p2 <- p1 +
  geom_fruit(geom = geom_tile, 
             mapping = aes(fill = Status, x = status_x), 
             pwidth = 0.2, 
             offset = 0.02,
             color = NA, 
             axis.params = list(
               axis = "x",
               text.size = 1.8,
               hjust = 0,
               vjust = 0.5,
               line.color = NA
             )) +
  scale_fill_manual(
    values = cols_status, 
    name = "LMM Status", 
    labels = c(
      "Enriched" = "Enriched in Sub-Saharan Africa",
      "Depleted" = "Depleted in Sub-Saharan Africa",
      "NS"       = "Not Significant"
    ),
    guide = guide_legend(
      direction = "horizontal", 
      title.position = "top", 
      order = 2
    )
  )

# Apply final theme 
p_final <- p2 +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",          
    legend.title = element_text(size = 6),          
    legend.text = element_text(size = 5),          
    legend.key.size = unit(2.5, "mm"),              
    legend.margin = margin(t = 2, b = 2),
    legend.box.margin = margin(t = 5, r = 0, b = 0, l = 0),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")
  )

# ==============================================================================
# 5. Save Output
# ==============================================================================

outfile_base <- "results/figures/Figure2c_PhylogeneticTree"

ggsave(paste0(outfile_base, ".pdf"), p_final, width = 130, height = 160, units = "mm", useDingbats = FALSE, bg = "transparent")
ggsave(paste0(outfile_base, ".png"), p_final, width = 130, height = 160, units = "mm", dpi = 300, bg = "white")

cat("Simplified tree plot successfully saved!\n")