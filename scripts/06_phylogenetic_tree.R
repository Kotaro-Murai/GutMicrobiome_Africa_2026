# ==============================================================================
# Script: 06_phylogenetic_tree.R
# Author: Kotaro Murai
# Description: 
#   Map LMM results onto the GTDB phylogenetic tree.
#   Features: Phylum ring, LMM Status ring, and 3 concentric rings showing 
#             mean abundance in Industrialized, Non-Ind, and Africa groups.
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

Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

# LMM Results from Script 05
res_df <- read_csv("results/tables/lmm_results_species_wide.csv", show_col_types = FALSE)

# Taxonomy data
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"
if("species" %in% colnames(taxa_data)) taxa_data <- taxa_data %>% rename(taxon_species = species)

# GTDB Tree & Metadata
tree_file <- "data/bac120_r207.tree"
meta_file <- "data/bac120_metadata_r207.tsv.gz" # Assuming compressed
if(!file.exists(tree_file)) stop("GTDB tree file not found in 'data/'!")

big_tree <- read.tree(tree_file)
gtdb_meta <- read_tsv(meta_file, col_select = c("accession", "gtdb_taxonomy", "gtdb_representative"), show_col_types = FALSE)

# Abundance & Metadata for concentric rings
species_data <- read_tsv("data/metalog_motus3_processed.tsv.gz", show_col_types = FALSE) %>% column_to_rownames(var = "rowname")
metadata <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)

# ==============================================================================
# 2. Data Matching & Tree Pruning
# ==============================================================================
cat("Matching mOTUs to GTDB tree and calculating mean abundances...\n")

# Clean GTDB metadata to get 's__Species_name' mapping
gtdb_map <- gtdb_meta %>%
  mutate(species_gtdb = str_extract(gtdb_taxonomy, "s__[^;]+$")) %>%
  arrange(species_gtdb, desc(gtdb_representative)) %>%
  distinct(species_gtdb, .keep_all = TRUE) %>%
  select(accession, species_gtdb)

# 1. Define LMM Status (Enriched / Depleted / NS)
res_clean <- res_df %>%
  filter(species != "-1") %>% # Remove unassigned
  mutate(
    motu_id = str_extract(species, "(?<=\\[).+?(?=\\]$)"),
    Status = case_when(
      q.value_AvsND < 0.05 & estimate_AvsND > 0.3 & q.value_AvsNDev < 0.05 & estimate_AvsNDev > 0.3 ~ "Enriched",
      q.value_AvsND < 0.05 & estimate_AvsND < -0.3 & q.value_AvsNDev < 0.05 & estimate_AvsNDev < -0.3 ~ "Depleted",
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
# 3. Calculate Mean Abundance per Group for Heatmap Rings
# ==============================================================================
cat("Calculating mean relative abundance per group...\n")

# Subset abundance data to only species present in the tree
tree_species_raw_names <- matched_data$species

# Calculate group means
group_abun <- metadata %>%
  select(sample_alias, Analysis_Group) %>%
  bind_cols(species_data[metadata$sample_alias, tree_species_raw_names]) %>%
  pivot_longer(cols = all_of(tree_species_raw_names), names_to = "species", values_to = "abundance") %>%
  group_by(Analysis_Group, species) %>%
  summarise(mean_abun = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(log10_abun = log10(mean_abun + 1E-5)) # Pseudo-count for log transformation

# Pivot wider for the 3 separate rings
abun_wide <- group_abun %>%
  pivot_wider(id_cols = species, names_from = Analysis_Group, values_from = log10_abun)

# Join everything into the final visualization dataframe
viz_data <- matched_data %>%
  select(accession, species, species_key, phylum, Status) %>%
  left_join(abun_wide, by = "species") %>%
  mutate(
    phylum_group = str_remove(phylum, "^p__"),
    phylum_group = replace_na(phylum_group, "Unassigned"),
    Status = factor(Status, levels = c("Enriched", "Depleted", "NS")),
    # Clean label for plotting (only for significant ones)
    clean_label = if_else(Status != "NS", 
                          str_replace_all(str_remove(species_key, "^s__"), "_", " "), 
                          NA_character_)
  )

d_tree <- data.frame(label = pruned_tree$tip.label) %>%
  left_join(viz_data, by = c("label" = "accession"))

# ==============================================================================
# 4. Plotting (Refined & Color-Synced)
# ==============================================================================
cat("Generating circular tree with concentric rings...\n")

heatmap_data <- viz_data %>%
  select(accession, Industrialized, `Non-Africa Non-Industrialized`, Africa) %>%
  pivot_longer(cols = -accession, names_to = "Group", values_to = "Abundance") %>%
  mutate(Group = factor(Group, levels = c("Industrialized", "Non-Africa Non-Industrialized", "Africa")))

p_base <- ggtree(pruned_tree, layout = "circular", size = 0.15, open.angle = 15)
max_radius <- max(p_base$data$x, na.rm = TRUE)
w_ring <- 0.08 * max_radius 

# --- Color Definitions ---
phylum_counts <- table(d_tree$phylum_group)
sorted_phyla <- names(sort(phylum_counts, decreasing = TRUE))
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

cols_status <- c(
  "Enriched" = thesis_colors[["Africa"]], 
  "Depleted" = thesis_colors[["Non-Africa Non-Industrialized"]],
  "NS"       = "#F5F5F5" 
)

# --- Build the Plot Layers ---
p <- p_base %<+% d_tree

# [Ring 1] Phylum
p1 <- p +
  geom_fruit(geom = geom_tile, mapping = aes(fill = phylum_group), pwidth = 0.08, offset = 0.04) +
  scale_fill_manual(values = cols_phylum, name = "Phylum", na.value = "gray95") +
  new_scale_fill()

# [Ring 2] Significance Status
p2 <- p1 +
  geom_fruit(geom = geom_tile, mapping = aes(fill = Status), pwidth = 0.08, offset = 0.04) +
  scale_fill_manual(values = cols_status, name = "LMM Status\n(Africa spec.)") +
  new_scale_fill()

# [Ring 3, 4, 5] Abundance Heatmaps
p3 <- p2 +
  geom_fruit(
    data = heatmap_data,
    geom = geom_tile,
    mapping = aes(y = accession, x = Group, fill = Abundance),
    pwidth = 0.12,
    offset = 0.06,
    color = NA
  ) +
  scale_fill_viridis_c(
    option = "magma", direction = -1, na.value = "white", 
    name = "Mean Abundance (Log10)\n[Inner: Ind. -> Mid: Non-Ind. -> Outer: Africa]"
  )

# [Ring 6] Text Labels
label_offset <- (w_ring * 5) + (max_radius * 0.25) # Push labels outside the rings

p_final <- p3 +
  geom_tiplab(
    aes(label = clean_label),
    size = 2.0,
    offset = label_offset,
    align = TRUE,
    linetype = "dotted",
    linesize = 0.2,
    color = "black"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(face = "bold", size = 10),
    legend.key.size = unit(0.5, "cm")
  )

# ==============================================================================
# 5. Save Output
# ==============================================================================
outfile_base <- "results/figures/Figure2b_PhylogeneticTree"

ggsave(paste0(outfile_base, ".pdf"), p_final, width = 22, height = 22, useDingbats = FALSE)
ggsave(paste0(outfile_base, ".png"), p_final, width = 22, height = 22, dpi = 300, bg = "white")

cat("Refined plot successfully saved!\n")