# ==============================================================================
# Script: 07_phylogenetic_tree_barplot.R
# Author: Kotaro Murai
# Description: 
#   Map LMM results onto the GTDB phylogenetic tree.
#   Features: Phylum ring, and Compressed Symmetrical Outward Barplot.
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

res_df <- read_csv("results/tables/lmm_results_species_wide.csv", show_col_types = FALSE)
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"
if("species" %in% colnames(taxa_data)) taxa_data <- taxa_data %>% rename(taxon_species = species)

tree_file <- "data/bac120_r207.tree"
meta_file <- "data/bac120_metadata_r207.tsv.gz"
if(!file.exists(tree_file)) stop("GTDB tree file not found in 'data/'!")

big_tree <- read.tree(tree_file)
gtdb_meta <- read_tsv(meta_file, col_select = c("accession", "gtdb_taxonomy", "gtdb_representative"), show_col_types = FALSE)

# ==============================================================================
# 2. Data Matching & Tree Pruning
# ==============================================================================
cat("Matching mOTUs to GTDB tree and extracting estimates...\n")

gtdb_map <- gtdb_meta %>%
  mutate(species_gtdb = str_extract(gtdb_taxonomy, "s__[^;]+$")) %>%
  arrange(species_gtdb, desc(gtdb_representative)) %>%
  distinct(species_gtdb, .keep_all = TRUE) %>%
  select(accession, species_gtdb)

res_clean <- res_df %>%
  filter(species != "-1") %>% 
  mutate(
    motu_id = str_extract(species, "(?<=\\[)[^\\[\\]]+(?=\\]$)"),
    Status_SSAvsInd = case_when(
      q.value_SSAvsInd < 0.05 & estimate_SSAvsInd > 0.1 ~ "Enriched",
      q.value_SSAvsInd < 0.05 & estimate_SSAvsInd < -0.1 ~ "Depleted",
      TRUE ~ "NS"
    ),
    Status_SSAvsONI = case_when(
      q.value_SSAvsONI < 0.05 & estimate_SSAvsONI > 0.1 ~ "Enriched",
      q.value_SSAvsONI < 0.05 & estimate_SSAvsONI < -0.1 ~ "Depleted",
      TRUE ~ "NS"
    )
  )

matched_step1 <- res_clean %>%
  inner_join(taxa_data, by = c("motu_id" = "id")) %>%
  mutate(species_key = if_else(str_detect(taxon_species, "^s__"), taxon_species, paste0("s__", taxon_species)))

matched_data <- matched_step1 %>%
  inner_join(gtdb_map, by = c("species_key" = "species_gtdb")) %>%
  filter(accession %in% big_tree$tip.label)

if(nrow(matched_data) == 0) stop("No species matched the tree.")
pruned_tree <- keep.tip(big_tree, matched_data$accession)

# ==============================================================================
# 3. Prepare Visualization Data 
# ==============================================================================
cat("Preparing tree metadata and calculating symmetric limits...\n")

viz_data <- matched_data %>%
  select(accession, species, species_key, phylum, 
         estimate_SSAvsInd, Status_SSAvsInd, 
         estimate_SSAvsONI, Status_SSAvsONI) %>%
  mutate(
    phylum_group = str_remove(phylum, "^p__"),
    phylum_group = replace_na(phylum_group, "Unassigned"),
    Status_SSAvsInd = factor(Status_SSAvsInd, levels = c("Enriched", "Depleted", "NS")),
    Status_SSAvsONI = factor(Status_SSAvsONI, levels = c("Enriched", "Depleted", "NS")),
    plot_est_Ind = if_else(Status_SSAvsInd == "NS", NA_real_, estimate_SSAvsInd),
    plot_est_ONI = if_else(Status_SSAvsONI == "NS", NA_real_, estimate_SSAvsONI)
  )

d_tree <- data.frame(label = pruned_tree$tip.label) %>%
  left_join(viz_data, by = c("label" = "accession"))

max_est <- max(abs(c(viz_data$plot_est_Ind, viz_data$plot_est_ONI)), na.rm = TRUE)
max_est <- ceiling(max_est * 10) / 10

dummy_idx <- which(d_tree$Status_SSAvsInd == "NS" & d_tree$Status_SSAvsONI == "NS")
if(length(dummy_idx) < 2) dummy_idx <- 1:2 

d_tree$plot_est_Ind[dummy_idx[1]] <- max_est
d_tree$plot_est_Ind[dummy_idx[2]] <- -max_est
d_tree$plot_est_ONI[dummy_idx[1]] <- 0.5
d_tree$plot_est_ONI[dummy_idx[2]] <- -0.5

d_tree$Status_SSAvsInd <- as.character(d_tree$Status_SSAvsInd)
d_tree$Status_SSAvsONI <- as.character(d_tree$Status_SSAvsONI)
d_tree$Status_SSAvsInd[dummy_idx[1:2]] <- "Dummy"
d_tree$Status_SSAvsONI[dummy_idx[1:2]] <- "Dummy"

d_tree$Status_SSAvsInd <- factor(d_tree$Status_SSAvsInd, levels = c("Enriched", "Depleted", "NS", "Dummy"))
d_tree$Status_SSAvsONI <- factor(d_tree$Status_SSAvsONI, levels = c("Enriched", "Depleted", "NS", "Dummy"))
# ----------------------------------------------------------------------

# ==============================================================================
# 4. Plotting (Phylum and LMM Estimate Barplot)
# ==============================================================================
cat("Generating circular tree with outward barplot...\n")

global_x_levels <- c("Phylum")
d_tree <- d_tree %>% mutate(phylum_x = factor("Phylum", levels = global_x_levels))

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
  "Dummy"    = "transparent" 
)

# --- Build the Plot Layers ---
p <- p_base %<+% d_tree

p1 <- p +
  geom_tippoint(
    mapping = aes(color = phylum_group),
    size = 0.5,  
    alpha = 0.9 
  ) +
  scale_color_manual(
    values = cols_phylum, 
    name = "Phylum", 
    na.value = "gray95",
    guide = guide_legend(
      direction = "horizontal", 
      nrow = 3, 
      title.position = "top", 
      order = 1,
      override.aes = list(size = 3) 
    ) 
  )

# [Ring 2 & 3] LMM Estimate Barplots
p2 <- p1 +
  geom_fruit(
    geom = geom_col,
    mapping = aes(x = plot_est_Ind, fill = Status_SSAvsInd),
    pwidth = 0.10,  
    offset = 0.07,  
    color = NA,
    orientation = "y",
    axis.params = list(
      axis = "x", text.size = 1.8, vjust = 0.5, hjust = -0.05,
      line.color = "black", line.size = 0.2,
      nbreak = 3
    ),
    grid.params = list(vline = FALSE, hline = FALSE, color = "grey75", size = 0.3)
  ) +
  geom_fruit(
    geom = geom_col,
    mapping = aes(x = plot_est_ONI, fill = Status_SSAvsONI),
    pwidth = 0.10,  
    offset = 0.05,
    color = NA,
    orientation = "y",
    axis.params = list(
      axis = "x", text.size = 1.8, vjust = 0.5, hjust = -0.05,
      line.color = "black", line.size = 0.2,
      nbreak = 3
    ),
    grid.params = list(vline = FALSE, hline = FALSE, color = "grey75", size = 0.3)
  ) +
  scale_fill_manual(
    values = cols_status, 
    name = "LMM Estimate", 
    labels = c(
      "Enriched" = "Enriched in Sub-Saharan Africa",
      "Depleted" = "Depleted in Sub-Saharan Africa"
    ),
    breaks = c("Enriched", "Depleted"), 
    na.translate = FALSE, 
    guide = guide_legend(direction = "horizontal", title.position = "top", order = 2)
  )

cat("Calculating MRCA nodes for target taxa highlights...\n")

tip_tax <- gtdb_meta %>% 
  filter(accession %in% pruned_tree$tip.label)
target_groups <- list(
  "Prevotella"              = "g__Prevotella\\b",
  "Gastranaerophilales"     = "o__Gastranaerophilales\\b",
  "Succinivibrio"           = "g__Succinivibrio\\b",
  "Bifidobacterium"         = "g__Bifidobacterium\\b",
  "Treponema"               = "g__Treponema_D\\b",
  "Phascolarctobacterium_A" = "g__Phascolarctobacterium_A\\b",
  "Bacteroides"             = "g__Bacteroides\\b",
  "Cryptobacteroides"       = "g__Cryptobacteroides\\b",
  "CAG_95"        = "g__CAG-95\\b",
  "CAG_83"        = "g__CAG-83\\b",
  "ER4"           = "g__ER4\\b",
  "Faecousia"     = "g__Faecousia\\b",
  "Acetatifactor" = "g__Acetatifactor\\b",
  "Enterocloster" = "g__Enterocloster\\b",
  "Lawsonibacter" = "g__Lawsonibacter\\b"
)
highlight_nodes <- purrr::map_dfr(names(target_groups), function(grp) {
  pattern <- target_groups[[grp]]
  matched_tips <- tip_tax %>% 
    filter(str_detect(gtdb_taxonomy, pattern)) %>% 
    pull(accession)
  if (length(matched_tips) == 0) return(NULL)
  node_id <- ifelse(
    length(matched_tips) == 1,
    which(pruned_tree$tip.label == matched_tips),
    MRCA(pruned_tree, matched_tips)
  )
  
  tibble(group = grp, node = node_id)
})

p3 <- p2 +
  new_scale_fill() + 
  geom_hilight(
    data = highlight_nodes,
    mapping = aes(node = node, fill = group),
    alpha = 0.3,
    extendto = 3,
    color = NA
  ) +
  scale_fill_manual(
    values = colorRampPalette(brewer.pal(8, "Set2"))(nrow(highlight_nodes)),
    name = "Highlighted Taxa"
  )

# Apply final theme 
p_final <- p3 +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",          
    legend.title = element_text(size = 6, face = "bold"),          
    legend.text = element_text(size = 5),          
    legend.key.size = unit(2.5, "mm"),              
    legend.margin = margin(t = 2, b = 2),
    legend.box.margin = margin(t = 5, r = 0, b = 0, l = 0),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")
  )

# ==============================================================================
# 5. Save Output
# ==============================================================================

outfile_base <- "results/figures/Figure2c_PhylogeneticTree_Barplot"

suppressWarnings({
  ggsave(paste0(outfile_base, ".pdf"), p_final, width = 140, height = 170, units = "mm", useDingbats = FALSE, bg = "transparent")
  ggsave(paste0(outfile_base, ".png"), p_final, width = 140, height = 170, units = "mm", dpi = 300, bg = "white")
})

cat("Compressed Barplot with grids successfully saved!\n")