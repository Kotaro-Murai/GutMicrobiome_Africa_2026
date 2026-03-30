# ==============================================================================
# Script: 00_plot_styles.R
# Description: Defines common themes and color palettes for the manuscript
# ==============================================================================

library(ggplot2)
library(ggthemes)
library(RColorBrewer)

# 2a. Master Color Palette for Phyla
phylum_master_colors <- c(
  "Firmicutes"        = "#FCCDE5",
  "Firmicutes_A"      = "#E8C5D0",
  "Firmicutes_B"      = "#D7C5D0",
  "Firmicutes_C"      = "#C28B99",
  "Bacteroidota"      = "#9DBAA6",
  "Proteobacteria"    = "#8C6BB1",
  "Actinobacteriota"  = "#E6C86E",
  "Verrucomicrobiota" = "#6BAED6",
  "Cyanobacteria"     = "#009E73",
  "Campylobacterota"  = "#FDB462",
  "Desulfobacterota"  = "#B3DE69",
  "Elusimicrobiota"   = "#BC80BD",
  "Synergistota"      = "#CCEBC5",
  "Fusobacteriota"    = "#FFED6F",
  "Spirochaetota"     = "#CC6677",
  "Others"            = "grey85",
  "Unassigned"        = "grey85"
)


# 2b. Color palette for microbiome composition (Phylum/Genus level bar plots)
pal_tableau <- tableau_color_pal("Tableau 20")
taxa_colors_20 <- pal_tableau(20)

# 3. Custom ggplot theme for the manuscript
theme_nature <- function(base_size = 6) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(size = base_size, color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size),
      plot.title = element_text(size = base_size + 1, hjust = 0), 
      
      line = element_line(linewidth = 0.25), 
      rect = element_rect(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      panel.border = element_rect(linewidth = 0.35, fill = NA), 
      
      panel.grid.major = element_line(linewidth = 0.15, color = "grey90"),
      panel.grid.minor = element_blank(),
      
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(3, "mm"), 
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
}


# 1. Main color palette for 3-group comparisons
group_colors <- c(
  "Industrialized" = "#56B4E9",        
  "Other Non-Industrialized" = "#009E73", 
  "Other\nNon-Industrialized" = "#009E73",
  "Sub-Saharan Africa" = "#E69F00"                             
)

# 4. Color palette for 2-group comparisons (Okabe-Ito style)
LMM_colors <- c(
  "Not significant" = "#9CA3AF",
  "Depleted"        = "#0072B2",
  "Enriched"        = "#D55E00"                                
)

# 5. Color palette for Lifestyle (4 groups)
lifestyle_colors <- c(
  "Hunter-gatherer"              = "#F5D44F", 
  "Rural_agrarian / traditional" = "#A3D362", 
  "Peri-urban / semi-urban"      = "#4CB140", 
  "Urban"                        = "#007A33"
)

# 5. Color palette for Lifestyle (4 groups)
lifestyle_colors <- c(
  "Hunter-gatherer"              = "#b3e2cd", 
  "Rural_agrarian / traditional" = "#fdcdac", 
  "Peri-urban / semi-urban"      = "#cbd5e8", 
  "Urban"                        = "#f4cae4"
)
