# ==============================================================================
# Script: 01_map_metagenome_data.R
# Author: Kotaro Murai
# Description: 
#   Visualize the global distribution of gut microbiome samples (Figure 1a).
#   Samples are aggregated by country, and circle size represents sample count.
# ==============================================================================

# Load libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
source("scripts/00_plot_styles.R")

# Create output directory if it doesn't exist
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Metadata and Aggregate ---

# Load processed metadata
# Ensure the file is placed in the 'data' directory
md <- read_tsv("data/metadata_with_groups_processed.tsv.gz") %>%
  mutate(
    Analysis_Group = factor(
      Analysis_Group, 
      levels = c("Sub-Saharan Africa", "Other Non-Industrialized", "Industrialized")
    )
  )

# Count samples per country
country_counts <- md %>%
  group_by(country, Analysis_Group) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

# --- 2. Prepare Map Data ---

# Load world map data (medium resolution)
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(name != "Antarctica")

# Standardize country names to match the map data
# (Mapping: Metadata Name -> rnaturalearth Name)
country_counts_fixed <- country_counts %>%
  mutate(country_mapped = case_when(
    country == "United States" ~ "United States of America",
    country == "Democratic Republic of the Congo" ~ "Dem. Rep. Congo",
    country == "Republic of Congo" ~ "Congo",
    country == "Cote d'Ivoire" ~ "Côte d'Ivoire",
    country == "Central African Republic" ~ "Central African Rep.",
    TRUE ~ country
  ))

# Merge map data with sample counts
map_data <- world %>%
  left_join(country_counts_fixed, by = c("name" = "country_mapped"))

# Calculate centroids for bubble placement
# Only calculate for countries with data to save processing time
map_data_centroids <- map_data %>%
  filter(!is.na(n)) %>%
  st_centroid()

# Manual correction for France
target_row <- which(map_data_centroids$name == "France")
if(length(target_row) > 0) {
  st_geometry(map_data_centroids)[[target_row]] <- st_point(c(2.2137, 46.2276))
}

# --- 3. Plotting (Figure 1a) ---

p_map <- ggplot() +
  # [Layer 1] Base World Map
  # Light gray land, white borders for a clean look
  geom_sf(data = world, 
          fill = "#E0E0E0", 
          color = "white", 
          linewidth = 0.05) +
  
  # [Layer 2] Sample Size Bubbles
  # Fixed color (Blue), size proportional to N
  geom_sf(data = map_data_centroids, 
          aes(size = n, fill = Analysis_Group), 
          shape = 21,        # Circle with border
          color = "white",   # Border color
          stroke = 0.1,      # Border width
          alpha = 0.8) +     # Transparency
  
  scale_fill_manual(values = group_colors, name = NULL, labels = function(x) stringr::str_wrap(x, width = 15)) +
  
  # Adjust bubble sizes
  scale_size_continuous(
    range = c(0.5, 5),           
    breaks = c(10, 100, 1000, 10000), 
    labels = scales::comma,    
    name = "Sample Size"
  ) +
  
  # Theme adjustments
  theme_void() + # Remove axes and background
  theme(
    text = element_text(family = "sans", size = 5, color = "black"), 
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(0, "mm"),
    legend.margin = margin(t = 0, b = 0),   
    legend.key.height = unit(1, "mm"),       
    legend.key.width = unit(1, "mm"), 
    legend.spacing.x = unit(0, "mm")
  ) + 
  guides(
    size = guide_legend(override.aes = list(fill = "grey50"))
  )

# Save plot
# Saving as PDF is recommended for vector quality in papers
ggsave("results/figures/Figure1a_Map_Overview.pdf", plot = p_map, 
       width = 65, height = 40, units = "mm", 
       useDingbats = FALSE)

# PNG
ggsave("results/figures/Figure1a_Map_Overview.png", plot = p_map, 
       width = 65, height = 40, units = "mm", dpi = 300)
