# ==============================================================================
# Script: 01_map_distribution.R
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
md <- read_tsv("data/metalog_metadata_processed.tsv.gz")

# Count samples per country
country_counts <- md %>%
  group_by(country) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

# --- 2. Prepare Map Data ---

# Load world map data (medium resolution)
world <- ne_countries(scale = "medium", returnclass = "sf")

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
# (Default centroid is shifted due to overseas territories like French Guiana)
target_row <- which(map_data_centroids$name == "France")
if(length(target_row) > 0) {
  st_geometry(map_data_centroids)[[target_row]] <- st_point(c(2.2137, 46.2276))
}

# --- 3. Plotting (Figure 1) ---

p_map <- ggplot() +
  # [Layer 1] Base World Map
  # Light gray land, white borders for a clean look
  geom_sf(data = world, 
          fill = "#F0F0F0", 
          color = "white", 
          size = 0.2) +
  
  # [Layer 2] Sample Size Bubbles
  # Fixed color (Blue), size proportional to N
  geom_sf(data = map_data_centroids, 
          aes(size = n), 
          shape = 21,        # Circle with border
          fill = "#377EB8",  # Uniform color
          color = "white",   # Border color
          stroke = 0.5,      # Border width
          alpha = 0.7) +     # Transparency
  
  # Adjust bubble sizes
  scale_size_area(
    max_size = 15,       # Adjust based on the actual figure size needs
    breaks = c(10, 100, 1000, 5000), 
    name = "Sample Size"
  ) +
  
  # Theme adjustments
  theme_void() + # Remove axes and background
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

# Display plot
print(p_map)

# Save plot
# Saving as PDF is recommended for vector quality in papers
ggsave("results/figures/Figure1a_Map_Overview.pdf", plot = p_map, width = 10, height = 6, useDingbats = FALSE)
ggsave("results/figures/Figure1a_Map_Overview.png", plot = p_map, width = 10, height = 6, dpi = 300)