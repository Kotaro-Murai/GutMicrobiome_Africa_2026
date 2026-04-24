# ==============================================================================
# Script:  Supple_plot_flagellar_heatmap.R
# Purpose: Identify and visualize flagellar features in Gastranaerophilales.
# Author:  Kotaro Murai
# ==============================================================================

# --- 0. Setup & Configuration ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(vegan)
  library(ggtext)
  library(KEGGREST)
})

# Load custom styles (ensure this file exists in your working directory)
source("scripts/00_plot_styles.R")

# Configuration List
CONFIG <- list(
  meta_file    = "data/genomes-all_metadata.tsv.gz",
  ko_matrix    = "data/genome_ko_matrix_species_rep_uhgg.csv.gz",
  mapping_file = "results/tables/flagellar_ko_gene_mapping.csv",
  out_pdf      = "results/figures/Supplementary_Figure3_Flagellar_Heatmap_Clustered.pdf",
  out_png      = "results/figures/Supplementary_Figure3_Flagellar_Heatmap_Clustered.png",
  target_order = "o__Gastranaerophilales",
  kegg_pathway = "path:map02040"
)

# Initialize directories
walk(c("results/figures", "results/tables"), ~dir.create(.x, showWarnings = FALSE, recursive = TRUE))
Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)


# ==============================================================================
# --- Utility Functions ---
# ==============================================================================

#' Fetch or Load Flagellar KO Mapping from KEGG
get_flagellar_mapping <- function(pathway_id, out_file) {
  if (file.exists(out_file)) {
    message("[-] Loading existing KEGG mapping table...")
    return(read_csv(out_file, show_col_types = FALSE))
  }
  
  message("[-] Fetching Flagellar KOs from KEGG...")
  # Get KO IDs for the pathway
  flag_kos <- sub("ko:", "", unname(keggLink("ko", pathway_id)))
  
  # Get all KO descriptions safely
  kegg_info <- tryCatch(
    { setNames(keggList("ko"), sub("ko:", "", names(keggList("ko")))) },
    error = function(e) { setNames(rep(NA, length(flag_kos)), flag_kos) }
  )
  
  # Parse KEGG strings (e.g., "fliC, hag; flagellin [KO:K02406]")
  parsed_info <- tibble(KO_ID = flag_kos) %>%
    mutate(
      raw_desc = kegg_info[KO_ID],
      Gene_Symbol = str_extract(raw_desc, "^[^;]+(?=;)") %>% str_split(",") %>% map_chr(~trimws(.x[1])),
      Description = str_extract(raw_desc, "(?<=;\\s)[^\\[]+") %>% str_remove_all("(?i)putative\\s+") %>% trimws()
    ) %>%
    replace_na(list(Gene_Symbol = "", Description = "Unknown")) %>%
    select(KO_ID, Gene_Symbol, Description)
  
  write_csv(parsed_info, out_file)
  message("[-] Saved mapping table to ", out_file)
  
  return(parsed_info)
}

#' Perform Hierarchical Clustering on Binary Matrix
cluster_elements <- function(mat, margin = 1, method_dist = "jaccard", method_hclust = "ward.D2") {
  target_mat <- if (margin == 1) mat else t(mat)
  dist_obj <- vegdist(target_mat, method = method_dist, binary = TRUE)
  hc_obj <- hclust(dist_obj, method = method_hclust)
  return(rownames(target_mat)[hc_obj$order])
}


# ==============================================================================
# --- Main Pipeline ---
# ==============================================================================

main <- function() {
  
  # 1. Load Data
  message("=== Step 1: Loading Data ===")
  
  # 1.1 Metadata
  uhgg_rep <- read_tsv(CONFIG$meta_file, show_col_types = FALSE) %>%
    filter(Genome == Species_rep) %>%
    mutate(
      domain  = str_extract(Lineage, "d__[^;]+"),
      order   = str_extract(Lineage, "o__[^;]+"),
      species = str_extract(Lineage, "s__[^;]+")
    ) %>%
    filter(domain == "d__Bacteria")
  
  target_genomes <- uhgg_rep %>%
    filter(order == CONFIG$target_order) %>%
    pull(Species_rep)
  
  # 1.2 KO Matrix
  message("[-] Loading KO matrix...")
  ko_matrix <- fread(CONFIG$ko_matrix) %>% as.data.frame()
  genome_col <- colnames(ko_matrix)[1]
  
  ko_mat_binary <- ko_matrix %>% 
    filter(!!sym(genome_col) %in% target_genomes) %>%
    column_to_rownames(var = genome_col) %>% 
    as.matrix() %>% 
    { +(. > 0) } # Convert to binary (0/1)
  
  # 2. Prepare Flagellar Mapping
  message("=== Step 2: Preparing Flagellar KO Mapping ===")
  ko_mapping <- get_flagellar_mapping(CONFIG$kegg_pathway, CONFIG$mapping_file)
  
  # Subset matrix to target KOs and remove empty columns
  valid_flag_kos <- intersect(ko_mapping$KO_ID, colnames(ko_mat_binary))
  flag_mat <- ko_mat_binary[, valid_flag_kos, drop = FALSE]
  flag_mat <- flag_mat[, colSums(flag_mat) > 0, drop = FALSE]
  
  # 3. Clustering & Plotting
  message("=== Step 3: Clustering and Plotting ===")
  
  # Order genomes (rows) and KOs (cols)
  genome_order <- cluster_elements(flag_mat, margin = 1)
  ko_order     <- cluster_elements(flag_mat, margin = 2)
  
  # Format data for ggplot
  heatmap_data <- as.data.frame(flag_mat) %>%
    rownames_to_column("Species_rep") %>%
    pivot_longer(cols = -Species_rep, names_to = "KO_ID", values_to = "Presence") %>%
    left_join(select(uhgg_rep, Species_rep, species), by = "Species_rep") %>%
    left_join(ko_mapping, by = "KO_ID") %>%
    mutate(
      Presence      = factor(Presence, levels = c(0, 1)),
      Species_Clean = str_replace(species, "s__", ""),
      Species_Label = sprintf("*%s* (%s)", Species_Clean, Species_rep),
      Species_rep   = factor(Species_rep, levels = rev(genome_order)),
      KO_ID         = factor(KO_ID, levels = ko_order)
    )
  
  # Generate Markdown Labels
  x_labels_map <- heatmap_data %>%
    distinct(KO_ID, Gene_Symbol) %>%
    mutate(Label = if_else(
      Gene_Symbol != "", 
      sprintf("*%s*<br><span style='font-size:4pt;'>(%s)</span>", Gene_Symbol, KO_ID), 
      KO_ID
    )) %>%
    pull(Label, name = KO_ID)
  
  y_labels_map <- heatmap_data %>%
    distinct(Species_rep, Species_Label) %>%
    pull(Species_Label, name = Species_rep)
  
  # Create Plot
  p_heatmap <- ggplot(heatmap_data, aes(x = KO_ID, y = Species_rep, fill = Presence)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_manual(
      values = c("0" = "grey95", "1" = "#4DBBD5"), 
      labels = c("Absent", "Present"), 
      name = "Gene Presence"
    ) +
    scale_y_discrete(labels = y_labels_map) +
    scale_x_discrete(labels = x_labels_map) +
    labs(
      title = "Presence/Absence of Flagellar KOs in Gastranaerophilales", 
      x = "Flagellar Genes (Clustered)", 
      y = "Species Representative"
    ) +
    theme_nature(base_size = 7) + 
    theme(
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      axis.text.y = element_markdown(size = 5),
      panel.grid  = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "bottom"
    )
  
  # Save Outputs
  ggsave(CONFIG$out_pdf, p_heatmap, width = 180, height = 220, units = "mm")
  ggsave(CONFIG$out_png, p_heatmap, width = 180, height = 220, units = "mm", dpi = 300)
  
  message("=== Pipeline Completed Successfully ===")
}

# Run the pipeline
main()