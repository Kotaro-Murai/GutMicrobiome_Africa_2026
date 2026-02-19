# ==============================================================================
# Script: 08_map_gastranaerophilales.R
# Description: Plot world map of Gastranaerophilales prevalence (Figure 3a).
# ==============================================================================

library(tidyverse)
library(data.table)
library(sf)
library(rnaturalearth)

source("scripts/00_plot_styles.R")
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---
cat("Loading data...\n")
Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

species_data <- read_tsv("data/metalog_motus3_processed.tsv.gz", show_col_types = FALSE) %>%
  column_to_rownames(var = "rowname")
metadata <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)
taxa_data <- read_tsv("data/mOTUs_3.0.0_GTDB_tax.tsv.gz", show_col_types = FALSE)

if(!"id" %in% colnames(taxa_data)) colnames(taxa_data)[1] <- "id"

# --- 2. Memory-Safe Aggregation Function---
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

# --- 3. Plot Function ---
plot_taxon_map <- function(target_rank, target_taxon_name, metric = "prevalence", min_n = 10, threshold = 0.0001) {
  
  cat(sprintf("Aggregating to '%s' level and calculating %s...\n", target_rank, metric))
  rank_data <- aggregate_by_rank(species_data, taxa_data, target_rank)
  
  if(!target_taxon_name %in% colnames(rank_data)) {
    stop(sprintf("Taxon '%s' not found in %s level data.", target_taxon_name, target_rank))
  }
  
  combined <- bind_cols(metadata, rank_data[, target_taxon_name, drop=FALSE]) %>%
    rename(abundance = all_of(target_taxon_name))
  
  country_stats <- combined %>%
    group_by(country) %>%
    summarise(
      n_samples = n(),
      value = if(metric == "abundance") {
        mean(abundance, na.rm = TRUE) * 100 
      } else {
        sum(abundance >= threshold, na.rm = TRUE) / n() * 100 
      },
      .groups = "drop"
    ) %>%
    filter(n_samples >= min_n)
  
  country_stats_fixed <- country_stats %>%
    mutate(country_map = case_when(
      country == "Central African Republic" ~ "Central African Rep.",
      country == "Cote d'Ivoire"  ~ "Côte d'Ivoire",
      country == "Democratic Republic of the Congo"  ~ "Dem. Rep. Congo",
      country == "Republic of Congo"  ~ "Congo",
      country == "United States"  ~ "United States of America",
      TRUE ~ country
    ))
  
  world_map <- ne_countries(scale = "medium", returnclass = "sf")
  map_data <- world_map %>% left_join(country_stats_fixed, by = c("name" = "country_map"))
  
  metric_label <- if(metric == "abundance") "Mean Abundance (%)" else "Prevalence (%)"
  display_name <- str_remove(target_taxon_name, "^[a-z]__")
  
  p <- ggplot(data = map_data) +
    geom_sf(aes(fill = value), color = "white", linewidth = 0.2) +
    
    scale_fill_distiller(
      palette = "YlOrRd",
      direction = 1,
      na.value = "grey90",
      name = metric_label,
      labels = function(x) paste0(x, "%")
    ) +
    
    labs(
      title = bquote(bold(.(ifelse(metric=="abundance", "Mean Abundance of", "Prevalence of"))) ~ bolditalic(.(display_name))),
      subtitle = paste0("Taxonomic Rank: ", str_to_title(target_rank), " | Samples >= ", min_n)
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40", margin = margin(b = 20)),
      legend.key.width = unit(2, "cm"),
      legend.title = element_text(vjust = 0.8, face = "bold")
    )
  
  return(p)
}

p_map <- plot_taxon_map("order", "o__Gastranaerophilales", metric = "prevalence", min_n = 10, threshold = 0.0001)

outfile <- "results/figures/Figure3a_Map_Gastranaerophilales"
ggsave(paste0(outfile, ".pdf"), p_map, width = 12, height = 7, useDingbats = FALSE)
ggsave(paste0(outfile, ".png"), p_map, width = 12, height = 7, dpi = 300, bg = "white")

cat("Successfully saved Figure 3a\n")