# ==============================================================================
# Script: 14_pcoa_cyanobacteria_uhgg.R
# Description: PCoA of bacterial genomes based on functional gene (KO) presence.
#              Highlights the Phylum Cyanobacteria (Figure 4a).
# ==============================================================================

# --- 0. Setup and Libraries ---
library(tidyverse)
library(data.table)
library(vegan)
source("scripts/00_plot_styles.R")

# Set directories
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Input Files
file_uhgg_meta <- "data/genomes-all_metadata.tsv.gz"
file_ko_matrix <- "data/genome_ko_matrix_species_rep_uhgg.csv.gz"

# ==============================================================================
# 1. UHGG Metadata & Taxonomy Preparation
# ==============================================================================
cat("=== Step 1: Preparing Metadata and Taxonomy ===\n")

Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

uhgg_data <- read_tsv(file_uhgg_meta, show_col_types = FALSE)

# Clean UHGG (Keep only Species Representatives)
uhgg_rep <- uhgg_data %>%
  filter(Genome == Species_rep) %>%
  mutate(
    gtdb_species = str_remove(str_extract(Lineage, "s__[^;]+"), "s__"),
    domain = str_extract(Lineage, "d__[^;]+"),
    phylum = str_extract(Lineage, "p__[^;]+")
  )

# ==============================================================================
# 2. Load Pre-computed KO Matrix
# ==============================================================================
cat("=== Step 2: Loading KO Matrix ===\n")

if (!file.exists(file_ko_matrix)) {
  stop(sprintf("Error: Pre-computed KO matrix not found at %s", file_ko_matrix))
}

ko_matrix <- data.table::fread(file_ko_matrix) %>% as.data.frame()

# ==============================================================================
# 3. PCoA Calculation (Jaccard Distance)
# ==============================================================================
cat("=== Step 3: Calculating PCoA ===\n")

# Filter for Bacteria only
bacterial_genomes <- uhgg_rep %>%
  filter(domain == "d__Bacteria") %>% 
  pull(Species_rep) %>%
  unique()

valid_genomes <- intersect(bacterial_genomes, ko_matrix$Genome)
cat(sprintf("Running ordination on %d bacterial genomes...\n", length(valid_genomes)))

# Prepare numeric matrix
bacteria_ko_numeric <- ko_matrix %>%
  filter(Genome %in% valid_genomes) %>%
  column_to_rownames("Genome") %>%
  as.matrix()

# Remove empty columns (KOs not present in this subset)
bacteria_ko_numeric <- bacteria_ko_numeric[, colSums(bacteria_ko_numeric) > 0]

# Calculate Jaccard distance (binary = TRUE ignores copy numbers, focusing on presence/absence)
dist_jaccard <- vegdist(bacteria_ko_numeric, method = "jaccard", binary = TRUE)

# Run PCoA
pcoa_res <- cmdscale(dist_jaccard, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Genome <- rownames(pcoa_df)

# Calculate explained variance
pcoa_eig <- pcoa_res$eig
pcoa_var <- pcoa_eig / sum(pcoa_eig)
pcoa_label_x <- sprintf("PCoA1 (%.1f%%)", pcoa_var[1] * 100)
pcoa_label_y <- sprintf("PCoA2 (%.1f%%)", pcoa_var[2] * 100)

# ==============================================================================
# 4. Visualization (Figure 4a)
# ==============================================================================
cat("=== Step 4: Generating Figure 4a ===\n")

phylum_counts <- uhgg_rep %>%
  filter(Species_rep %in% valid_genomes) %>%
  mutate(Phylum_clean = str_remove(phylum, "p__")) %>%
  count(Phylum_clean, sort = TRUE)

top_n <- 10
top_phyla <- head(phylum_counts$Phylum_clean, top_n)

# Combine PCoA coordinates with Taxonomy and prepare highlighting
plot_data <- pcoa_df %>%
  inner_join(uhgg_rep %>% select(Genome = Species_rep, phylum), by = "Genome") %>%
  mutate(
    Phylum_clean = str_remove(phylum, "p__"),
    Phylum_Group = if_else(Phylum_clean %in% top_phyla, Phylum_clean, "Other Phyla"),
    Phylum_Group = factor(Phylum_Group, levels = c(sort(top_phyla), "Other Phyla"))
  ) %>%
  mutate(Sort_Order = case_when(
    Phylum_Group == "Other Phyla" ~ 1,
    Phylum_Group == "Cyanobacteria" ~ 3,
    TRUE ~ 2
  )) %>%
  arrange(Sort_Order)

existing_master <- phylum_master_colors[names(phylum_master_colors) %in% top_phyla]
missing_names <- setdiff(top_phyla, names(phylum_master_colors))

if(length(missing_names) > 0) {
  extra_cols <- taxa_colors_20[1:length(missing_names)]
  names(extra_cols) <- missing_names
  final_palette <- c("Other Phyla" = "grey85", existing_master, extra_cols)
} else {
  final_palette <- c("Other Phyla" = "grey85", existing_master)
}

# Plot Figure 4a
p_fig4a <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2, color = Phylum_Group)) +
  
  geom_vline(xintercept = 0, color = "gray90", linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = 0, color = "gray90", linetype = "dashed", linewidth = 0.25) +
  stat_ellipse(
    data = filter(plot_data, Phylum_Group == "Cyanobacteria"),
    aes(group = Phylum_Group, color = Phylum_Group),
    geom = "polygon", 
    level = 0.99,     
    fill = NA,        
    linewidth = 0.3,   
    linetype = "solid", 
    alpha = 0.1, 
    show.legend = FALSE
  ) +
  geom_point(aes(alpha = ifelse(Phylum_Group == "Other Phyla", 0.5, 0.5)), 
             size = 0.8, stroke = 0) +
  
  scale_color_manual(values = final_palette) +
  scale_alpha_identity() +
  
  labs(
    x = pcoa_label_x,
    y = pcoa_label_y,
    color = "Phylum"
  ) +
  
  theme_nature() +
  theme(
    legend.position = c(0.99, 0.01), 
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.key.size = unit(2, "mm"),
    legend.direction = "vertical",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6, face = "bold"),
    legend.margin = margin(l = 1),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
  ) +
  guides(color = guide_legend(ncol = 1, 
                              byrow = TRUE,
                              title.position = "top",
                              override.aes = list(size = 1.5, alpha = 1)
                              ))

# Save outputs
outfile <- "results/figures/Figure4a_PCoA_UHGG_Cyanobacteria"

ggsave(paste0(outfile, ".pdf"), p_fig4a, width = 90, height = 82, units = "mm", useDingbats = FALSE)
ggsave(paste0(outfile, ".png"), p_fig4a, width = 90, height = 82, units = "mm", dpi = 300, bg = "white")

cat("Successfully saved Figure 4a!\n")

# ==============================================================================
# 6. Export Supplementary Table 12 (Bacterial Genomes Information)
# ==============================================================================
cat("=== Step 6: Exporting Supplementary Table 12 (All Valid Genomes) ===\n")

supp_table_12 <- uhgg_rep %>%
  filter(Species_rep %in% valid_genomes) %>%
  mutate(
    Domain  = str_remove(str_extract(Lineage, "d__[^;]+"), "d__"),
    Phylum  = str_remove(str_extract(Lineage, "p__[^;]+"), "p__"),
    Class   = str_remove(str_extract(Lineage, "c__[^;]+"), "c__"),
    Order   = str_remove(str_extract(Lineage, "o__[^;]+"), "o__"),
    Family  = str_remove(str_extract(Lineage, "f__[^;]+"), "f__"),
    Genus   = str_remove(str_extract(Lineage, "g__[^;]+"), "g__"),
    Species = str_remove(str_extract(Lineage, "s__[^;]+"), "s__")
  ) %>%
  select(
    Species_Representative = Species_rep,
    Genome_Type = Genome_type,
    Domain,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Length,
    N_contigs,
    N50,
    GC_content,
    Completeness,
    Contamination,
    Country,
    Continent,
    Study_Accession = Study_accession,
    Sample_Accession = Sample_accession
  ) %>%
  arrange(Phylum, Class, Order, Family, Genus, Species)

write_csv(supp_table_12, "results/tables/Supplementary_Table12_Bacterial_Genomes_UHGG.csv")

cat("Successfully saved Supplementary Table 12 (Bacterial Genomes Info)!\n")