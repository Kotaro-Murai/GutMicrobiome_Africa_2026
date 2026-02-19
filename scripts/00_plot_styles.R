# ==============================================================================
# Script: 00_plot_styles.R
# Description: Defines common themes and color palettes for the manuscript
# ==============================================================================

library(ggplot2)
library(ggthemes)
library(RColorBrewer)

# 1. Main color palette for 3-group comparisons
thesis_colors <- c(
  "Industrialized" = "#9CA3AF",      
  "Non-Africa Non-Industrialized" = "#2874A6", 
  "Africa" = "#E67E22"                            
)

# 2a. Master Color Palette for Phyla
phylum_master_colors <- c(
  "Firmicutes"        = "#D9A5B3",
  "Firmicutes_A"      = "#E8C5D0",
  "Firmicutes_C"      = "#C28B99",
  "Bacteroidota"      = "#9DBAA6",
  "Proteobacteria"    = "#8C6BB1",
  "Actinobacteriota"  = "#E6C86E",
  "Verrucomicrobiota" = "#6BAED6",
  "Cyanobacteria"     = "#009E73",
  "Campylobacterota"  = "#FDB462",
  "Desulfobacterota"  = "#B3DE69",
  "Spirochaetota"     = "#FCCDE5",
  "Elusimicrobiota"   = "#BC80BD",
  "Synergistota"      = "#CCEBC5",
  "Fusobacteriota"    = "#FFED6F",
  "Others"            = "grey85",
  "Unassigned"        = "grey85"
)


# 2b. Color palette for microbiome composition (Phylum/Genus level bar plots)
pal_tableau <- tableau_color_pal("Tableau 20")
taxa_colors_20 <- pal_tableau(20)

# 3. Custom color palette (for alternative use cases)
tidy_col <- c("#18BC9C","#CCBE93","#a6cee3","#1f78b4","#b2df8a",
              "#fb9a99","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")

# 4. Custom ggplot theme for the manuscript
theme_thesis <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      axis.text = element_text(color = "black")
    )
}