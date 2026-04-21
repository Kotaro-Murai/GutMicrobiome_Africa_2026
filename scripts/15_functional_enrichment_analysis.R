# ==============================================================================
# Script: 15_functional_enrichment_gastranaerophilales.R
# Description: Identify Gastranaerophilales-specific functional features.
#              Figure 4b: KEGG Pathway Enrichment Analysis
#              Figure 4c: Top Enriched KOs (Prevalence Comparison)
# ==============================================================================

# --- 0. Setup and Libraries ---
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(KEGGREST)
library(viridis)
library(ggrepel)
library(patchwork)
library(ggtext)

# Load custom styles
source("scripts/00_plot_styles.R")

# Set directories
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# Input Files
file_uhgg_meta <- "data/genomes-all_metadata.tsv.gz"
file_ko_matrix <- "data/genome_ko_matrix_species_rep_uhgg.csv.gz"

# ==============================================================================
# 1. Load Data & Define Groups
# ==============================================================================
cat("=== Step 1: Loading Data and Defining Groups ===\n")

Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

# Load UHGG metadata (Add Order extraction)
uhgg_data <- read_tsv(file_uhgg_meta, show_col_types = FALSE)
uhgg_rep <- uhgg_data %>%
  filter(Genome == Species_rep) %>%
  mutate(
    domain = str_extract(Lineage, "d__[^;]+"),
    order  = str_extract(Lineage, "o__[^;]+")
  ) %>%
  filter(domain == "d__Bacteria")

valid_genomes <- uhgg_rep$Species_rep

# Load KO Matrix
cat("Loading KO matrix...\n")
ko_matrix <- data.table::fread(file_ko_matrix) %>% as.data.frame()
genome_col <- colnames(ko_matrix)[1]

# Filter matrix to ONLY include genomes present in this study
ko_matrix <- ko_matrix %>% filter(!!sym(genome_col) %in% valid_genomes)

# Define Groups
target_genomes <- uhgg_rep %>%
  filter(order == "o__Gastranaerophilales") %>%
  pull(Species_rep)

group_vec <- if_else(ko_matrix[[genome_col]] %in% target_genomes, "Gastranaerophilales", "Other Bacteria")

# Subset Matrix for Analysis
ko_mat <- ko_matrix %>% column_to_rownames(var = genome_col) %>% as.matrix()
rownames(ko_mat) <- ko_matrix[[genome_col]]

# Remove KOs completely absent in this cohort
ko_mat <- ko_mat[, colSums(ko_mat > 0) > 0]

cat(sprintf("Cohort restricted. Comparing %d Gastranaerophilales genomes vs %d Other genomes.\n", 
            sum(group_vec == "Gastranaerophilales"), sum(group_vec == "Other Bacteria")))

# ==============================================================================
# 2. Statistical Analysis (Wilcoxon & Prevalence Difference)
# ==============================================================================
cat("=== Step 2: Statistical Testing ===\n")

stats_df <- data.frame(row.names = colnames(ko_mat))
ko_mat_binary <- +(ko_mat > 0)

# A) Prevalence (Mean Abundance in binary matrix)
mean_gastra <- colMeans(ko_mat_binary[group_vec == "Gastranaerophilales", , drop=FALSE])
mean_other  <- colMeans(ko_mat_binary[group_vec == "Other Bacteria", , drop=FALSE])

stats_df$Prevalence_Gastra <- mean_gastra
stats_df$Prevalence_Other  <- mean_other
stats_df$Diff_Prev <- mean_gastra - mean_other

# B) Wilcoxon Rank Sum Test
ko_to_test <- names(mean_gastra)[mean_gastra > 0.1]
cat(sprintf("Testing %d KOs (Prevalence > 10%% in Gastranaerophilales)...\n", length(ko_to_test)))

p_vals <- apply(ko_mat_binary[, ko_to_test, drop=FALSE], 2, function(x) {
  if(var(x) == 0) return(1)
  wilcox.test(x ~ group_vec, alternative = "greater", 
              subset = group_vec %in% c("Gastranaerophilales", "Other Bacteria"))$p.value
})

stats_df$p_value <- NA
stats_df[names(p_vals), "p_value"] <- p_vals
stats_df$q_value <- p.adjust(stats_df$p_value, method = "BH")


# ==============================================================================
# Filter Significant KOs & Export Supplementary Table 14
# ==============================================================================
sig_threshold_q <- 0.05
sig_threshold_diff <- 0.3

target_kos <- rownames(stats_df)[which(stats_df$q_value < sig_threshold_q & stats_df$Diff_Prev > sig_threshold_diff)]
background_kos <- rownames(stats_df)

cat(sprintf(" -> Found %d significant Gastranaerophilales-enriched KOs.\n", length(target_kos)))

cat("=== Exporting Supplementary Table 14 (All KO Stats) ===\n")

supp_table_14 <- stats_df %>% 
  rownames_to_column("KO_ID") %>%
  select(
    KO_ID,
    Prevalence_Gastranaerophilales = Prevalence_Gastra,
    Prevalence_Other_Bacteria = Prevalence_Other,
    Prevalence_Difference = Diff_Prev,
    P_value = p_value,
    FDR = q_value
  ) %>%
  arrange(FDR, desc(Prevalence_Difference))

write_csv(supp_table_14, "results/tables/Supplementary_Table14_All_KOs_Stats.csv")
cat("Successfully saved Supplementary Table 14 (All KO Stats)!\n")

# ==============================================================================
# 3. Figure 4b: KEGG Pathway Enrichment Analysis & Supp Table 13
# ==============================================================================
cat("=== Step 3: KEGG Enrichment (Figure 4b & Supp Table 13) ===\n")

if(length(target_kos) > 0) {
  
  tryCatch({
    link_ko_path <- keggLink("pathway", "ko")
    map_ko_path <- data.frame(pathway_id = sub("path:", "", unname(link_ko_path)), 
                              ko_id = sub("ko:", "", names(link_ko_path)))
    
    list_path <- keggList("pathway")
    map_path_name <- data.frame(pathway_id = sub("path:", "", names(list_path)), 
                                pathway_name = unname(list_path))
    
    max_kegg_size <- 500 
    
    pathway_sizes <- map_ko_path %>%
      filter(str_detect(pathway_id, "^map")) %>%
      count(pathway_id, name = "kegg_size")
    
    valid_pathways <- pathway_sizes %>%
      filter(kegg_size <= max_kegg_size) %>%
      pull(pathway_id)
    
    map_ko_path_clean <- map_ko_path %>% 
      filter(pathway_id %in% valid_pathways)
    
    map_path_name_clean <- map_path_name %>% 
      filter(pathway_id %in% valid_pathways)
    # --------------------------------------------------------------------------
    ek <- enricher(
      gene = target_kos,
      universe = background_kos,
      TERM2GENE = map_ko_path_clean,
      TERM2NAME = map_path_name_clean,
      pAdjustMethod = "BH",
      pvalueCutoff = 1.0,  # All p-values
      qvalueCutoff = 1.0,  # All q-values
      minGSSize = 5,
      maxGSSize = 500 
    )
    
    if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
      res_df <- as.data.frame(ek) %>%
        separate(GeneRatio, c("k", "n"), sep = "/", remove = FALSE) %>%
        separate(BgRatio, c("M", "N"), sep = "/", remove = FALSE) %>%
        mutate(
          k = as.numeric(k), n = as.numeric(n),
          M = as.numeric(M), N = as.numeric(N),
          FoldEnrichment = (k / n) / (M / N)
        ) %>%
        select(-k, -n, -M, -N) %>% 
        arrange(pvalue) 
      write_csv(res_df, "results/tables/Supplementary_Table13_KEGG_Enrichment_All.csv")
      cat("Successfully saved Supplementary Table 13 (All KEGG Results)!\n")
      
      plot_data_b <- res_df %>%
        filter(p.adjust < 0.05) %>%
        mutate(LogQ = -log10(p.adjust)) %>%
        arrange(FoldEnrichment) %>%
        slice_tail(n = 10) %>% 
        mutate(Description = factor(Description, levels = unique(Description)))
      
      if(nrow(plot_data_b) > 0) {
        p_fig4b <- ggplot(plot_data_b, aes(x = FoldEnrichment, y = Description)) +
          geom_vline(xintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.2) +
          geom_point(aes(size = Count, color = LogQ)) +
          scale_color_distiller(palette = "YlOrRd", direction = 1, name = expression("-Log"[10]*"(FDR)")) +
          scale_size_continuous(range = c(1, 3), name = "Count") +
          scale_x_continuous(limits = c(min(0.8, min(plot_data_b$FoldEnrichment)), NA), 
                             expand = expansion(mult = c(0, 0.1))) +
          labs(x = "Fold Enrichment", y = NULL) + 
          theme_nature(base_size = 6) +
          theme(
            axis.text.y = element_text(size = 5),
            legend.position = "bottom", 
            legend.key.height = unit(2, "mm"),
            legend.text = element_text(size = 4),
            legend.title = element_text(size = 5, face = "bold")
          )
      } else {
        cat("No significant pathways (FDR < 0.05) found for plotting.\n")
      }
      
    } else {
      cat("No pathways found.\n")
    }
  }, error = function(e) {
    cat("Error in KEGG enrichment: ", e$message, "\n")
  })
}

# ==============================================================================
# 4. Figure 4c: Top Enriched KOs
# ==============================================================================
cat("=== Step 4: Top Enriched KOs (Figure 4c) ===\n")

top_kos_df <- stats_df %>%
  rownames_to_column("KO") %>%
  filter(q_value < 0.05) %>% 
  arrange(desc(Diff_Prev)) %>%
  slice_head(n = 10)

if(nrow(top_kos_df) > 0) {
  
  manual_label_file <- "data/manual_ko_labels.csv"
  
  if(file.exists(manual_label_file)) {
    cat("Reading manually edited KO labels...\n")
    manual_labels <- read_csv(manual_label_file, show_col_types = FALSE)
    top_kos_df <- top_kos_df %>%
      left_join(manual_labels %>% select(KO, Manual_Label), by = "KO") %>%
      mutate(Label = Manual_Label)
  } else {
    cat("Fetching KO descriptions from KEGG...\n")
    ko_list <- top_kos_df$KO
    ko_desc <- tryCatch({
      kegg_info <- keggList("ko")
      names(kegg_info) <- sub("ko:", "", names(kegg_info))
      kegg_info[ko_list]
    }, error = function(e) {
      cat("KEGG API Error. Using IDs only.\n")
      return(setNames(ko_list, ko_list))
    })
    clean_desc <- sapply(ko_desc, function(x) {
      if (is.na(x)) return("Unknown")
      
      parts <- str_split(x, ";", n = 2)[[1]]
      
      if (length(parts) == 2) {
        gene_symbol <- trimws(parts[1])
        desc_part <- str_split(trimws(parts[2]), "\\[")[[1]][1]
      } else {
        gene_symbol <- ""
        desc_part <- str_split(x, "\\[")[[1]][1]
      }
      desc_part <- str_replace_all(desc_part, "(?i)putative\\s+", "")
      desc_part <- trimws(desc_part)
      if (nchar(gene_symbol) > 0) {
        full_label <- paste0("*", gene_symbol, "* (", desc_part, ")")
      } else {
        full_label <- desc_part
      }
      wrapped_label <- str_wrap(full_label, width = 70)
      html_label <- str_replace_all(wrapped_label, "\n", "<br>")
      
      return(html_label)
    })
    
    top_kos_df$Label <- clean_desc[top_kos_df$KO]
    
    top_kos_df %>% 
      select(KO, Label) %>% 
      rename(Manual_Label = Label) %>%
      write_csv("data/manual_ko_labels_template.csv")
    cat("Saved template for manual labels to 'data/manual_ko_labels_template.csv'.\n")
  }
  
  plot_data_c <- top_kos_df %>%
    select(Label, Prevalence_Gastra, Prevalence_Other) %>%
    pivot_longer(cols = c(Prevalence_Gastra, Prevalence_Other), 
                 names_to = "Group", values_to = "Prevalence") %>%
    mutate(
      Group = if_else(Group == "Prevalence_Gastra", "Gastranaerophilales", "Other Bacteria"),
      Label = factor(Label, levels = rev(unique(top_kos_df$Label))) 
    )
  
  p_fig4c <- ggplot(plot_data_c, aes(x = Prevalence * 100, y = Label, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             width = 0.7, color = "black", linewidth = 0.1) +
    scale_fill_manual(values = c("Gastranaerophilales" = "#008B8B", "Other Bacteria" = "grey80")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)), breaks = c(0, 50, 100)) +
    labs(x = "Prevalence (%)", y = NULL) + 
    theme_nature(base_size = 6) +
    theme(
      axis.text.y = element_markdown(size = 5),
      legend.position = "bottom", 
      legend.title = element_blank(),
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 5)
    )
  
  if(exists("p_fig4b")) {
    cat("Combining Figure 4b and 4c...\n")
    combined_plot <- (p_fig4b / p_fig4c) + 
      plot_layout(heights = c(1, 1)) 
    
    ggsave("results/figures/Figure4bc_Combined.pdf", combined_plot, width = 90, height = 90, units = "mm")
    ggsave("results/figures/Figure4bc_Combined.png", combined_plot, width = 90, height = 90, units = "mm", dpi = 300)
    cat("Successfully saved combined Figure 4bc!\n")
    
  } else {
    ggsave("results/figures/Figure4c_Top_KOs_Barplot_Gastranaerophilales.pdf", p_fig4c, width = 90, height = 65, units = "mm")
  }
  
} else {
  cat("No highly specific KOs found (Diff > 0.3).\n")
}