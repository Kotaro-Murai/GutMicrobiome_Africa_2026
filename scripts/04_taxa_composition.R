# ==============================================================================
# Script: 04_taxa_composition.R
# ==============================================================================

# Load libraries
library(tidyverse)
library(data.table)
library(patchwork)
library(ggh4x)
source("scripts/00_plot_styles.R")

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---
cat("Loading data...\n")
Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

species_data <- read_tsv("data/metalog_motus3_processed.tsv.gz") %>% column_to_rownames(var = "rowname")
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)
if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"
metadata <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)

if("-1" %in% colnames(species_data)) {
  species_data <- species_data %>% select(-`-1`)
}

metadata <- metadata %>%
  mutate(Analysis_Group = case_when(
    Analysis_Group == "Other Non-Industrialized" ~ "Other Non-Industrialized",
    TRUE ~ Analysis_Group
  ))

metadata$Analysis_Group <- factor(metadata$Analysis_Group, 
                                  levels = c("Sub-Saharan Africa", "Other Non-Industrialized", "Industrialized"))

strip_colors <- scales::alpha(unname(group_colors[levels(metadata$Analysis_Group)]), 0.7)

ranks_to_plot <- c("phylum", "genus")
top_n_taxa <- 10
plot_list <- list()

for (target_rank in ranks_to_plot) {
  
  cat("Processing rank:", target_rank, "...\n")
  
  col_names <- colnames(species_data)
  ids <- str_extract(col_names, "(?<=\\[)[^\\[\\]]+(?=\\]$)")
  
  map_df <- data.frame(original_col = col_names, id = ids) %>%
    left_join(taxa_data, by = "id") %>%
    mutate(
      RankName = .[[target_rank]],
      RankName = if_else(is.na(RankName) | RankName == "", "Unassigned", RankName),
      CleanName = str_remove(RankName, "^[pcofg]__") 
    ) %>%
    filter(!is.na(id))
  
  valid_cols <- map_df$original_col
  country_species_df <- bind_cols(metadata %>% select(country, Analysis_Group), species_data[, valid_cols]) %>%
    mutate(
      country = case_when(
        country == "Democratic Republic of the Congo" ~ "Dem. Rep. Congo",
        country == "Republic of Congo" ~ "Rep. Congo",
        country == "Central African Republic" ~ "Central African Rep.",
        country == "United States" ~ "USA",
        country == "United Kingdom" ~ "UK",
        TRUE ~ country
      )
    ) %>%
    group_by(country, Analysis_Group) %>%
    summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  plot_source <- country_species_df %>%
    pivot_longer(cols = -c(country, Analysis_Group), names_to = "original_col", values_to = "MeanAbundance") %>%
    left_join(map_df %>% select(original_col, CleanName), by = "original_col") %>%
    group_by(country, Analysis_Group, CleanName) %>%
    summarise(MeanAbundance = sum(MeanAbundance, na.rm = TRUE), .groups = "drop") %>%
    rename(Taxon = CleanName) %>%
    group_by(country) %>%
    mutate(RelAbundance = MeanAbundance / sum(MeanAbundance)) %>%
    ungroup()
  
  # --- Top Taxa ---
  top_taxa_list <- plot_source %>%
    filter(Taxon != "Unassigned") %>%
    group_by(Taxon) %>%
    summarise(TotalMean = mean(RelAbundance)) %>%
    arrange(desc(TotalMean)) %>%
    slice_head(n = top_n_taxa) %>%
    pull(Taxon)
  
  custom_levels <- unique(c("Bacteroides", "Prevotella", top_taxa_list))
  custom_levels <- custom_levels[custom_levels %in% top_taxa_list]
  custom_levels <- c(custom_levels, "Others") 
  
  plot_data <- plot_source %>%
    mutate(
      DisplayTaxon = if_else(Taxon %in% top_taxa_list, Taxon, "Others"),
      DisplayTaxon = factor(DisplayTaxon, levels = custom_levels)
    ) %>%
    group_by(country, Analysis_Group, DisplayTaxon) %>%
    summarise(RelAbundance = sum(RelAbundance), .groups = "drop")
  
  # --- Sorting ---
  group_dominant_taxa <- plot_data %>%
    filter(DisplayTaxon != "Others") %>%
    group_by(Analysis_Group, DisplayTaxon) %>%
    summarise(MeanAb = mean(RelAbundance), .groups = "drop") %>%
    group_by(Analysis_Group) %>%
    slice_max(MeanAb, n = 1) %>% 
    select(Analysis_Group, DominantTaxon = DisplayTaxon)
  
  country_order <- plot_data %>%
    left_join(group_dominant_taxa, by = "Analysis_Group") %>%
    filter(DisplayTaxon == DominantTaxon) %>% 
    arrange(Analysis_Group, desc(RelAbundance)) %>% 
    pull(country)
  
  plot_data$country <- factor(plot_data$country, levels = country_order)
  
  plot_taxa <- levels(plot_data$DisplayTaxon)
  my_colors <- character(length(plot_taxa))
  names(my_colors) <- plot_taxa
  fallback_idx <- 1 
  
  for (tax in plot_taxa) {
    if (tax %in% names(phylum_master_colors)) {
      my_colors[tax] <- phylum_master_colors[tax]
    } else if (tax == "Others") {
      my_colors[tax] <- "grey90"
    } else {
      my_colors[tax] <- taxa_colors_20[fallback_idx]
      fallback_idx <- fallback_idx + 1
      if(fallback_idx > length(taxa_colors_20)) fallback_idx <- 1
    }
  }
  
  p_base <- ggplot(plot_data, aes(x = country, y = RelAbundance, fill = DisplayTaxon)) +
    geom_bar(stat = "identity", width = 0.95, color = NA) + 
    scale_fill_manual(values = my_colors, name = str_to_title(target_rank)) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(title = NULL, y = "Relative Abundance", x = NULL) +
    theme_nature(base_size = 6) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.key.size = unit(2.5, "mm"),    
      legend.text = element_text(size = 5), 
      legend.title = element_text(size = 6),
      legend.margin = margin(l = 2)
    ) +
    guides(fill = guide_legend(ncol = 1))
  
  if (target_rank == "phylum") {
    p <- p_base + 
      facet_grid2(~ Analysis_Group, scales = "free_x", space = "free_x",
                  strip = strip_themed(background_x = elem_list_rect(fill = strip_colors))) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 6, color = "white", margin = margin(2, 0, 2, 0)),
        strip.background = element_rect(color = "black", linewidth = 0.25)
      )
  } else {
    p <- p_base + 
      facet_grid(~ Analysis_Group, scales = "free_x", space = "free_x") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5, color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.justification = c(0, 0),          
        legend.box.margin = margin(b = -30)
      )
  }
  
  plot_list[[target_rank]] <- p
}

cat("\nCombining and saving plots...\n")

p_combined <- plot_list$phylum / plot_list$genus + 
  plot_layout(heights = c(1, 1))

filename_base <- "results/figures/Figure1d_Taxa_Composition"

ggsave(paste0(filename_base, ".pdf"), p_combined, width = 180, height = 75, units = "mm", useDingbats = FALSE, bg = "transparent")
ggsave(paste0(filename_base, ".png"), p_combined, width = 180, height = 75, units = "mm", dpi = 300, bg = "white")

cat(sprintf("Successfully saved stacked plots to %s(.pdf/.png)\n", filename_base))