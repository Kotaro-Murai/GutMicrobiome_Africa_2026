# ==============================================================================
# Script: 02_pca_country.R
# Author: Kotaro Murai
# Description: 
#   Perform Principal Component Analysis (PCA) on species-level relative abundances.
#   Define three analysis groups (Industrialized, Other Non-Industrialized, 
#   and Sub-Saharan Africa) based on the PC2 score threshold (Malaysia).
# ==============================================================================

# Load libraries
library(tidyverse)
library(ggrepel)

source("scripts/00_plot_styles.R")

# Create output directories if they don't exist
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---

Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

# Read compressed abundance table and metadata
d <- read_tsv("data/metalog_motus3_processed.tsv.gz") %>%
  column_to_rownames(var = "rowname")
md <- read_tsv("data/metalog_metadata_processed.tsv.gz")

# --- 2. Filtering for Analysis ---

# Filter species: mean relative abundance > 0.01% (1E-4) AND prevalence > 10% (0.1)
keep <- apply(d, 2, mean) > 1E-4 & apply(d > 0, 2, mean) > 0.1
d_filtered <- d[, keep]

# --- 3. Principal Component Analysis (PCA) ---

# Log10 transformation with pseudo-count (1E-4)
d.pr <- prcomp(log10(d_filtered + 1E-4))

# Calculate explained variance for axis labels
ex_var <- summary(d.pr)$importance[2, ]
x.imp <- round(ex_var[1], digits = 3) * 100
y.imp <- round(ex_var[2], digits = 3) * 100
xlab <- paste0("PC1 (", x.imp, "%)")
ylab <- paste0("PC2 (", y.imp, "%)")

# --- 4. Define Groups based on Malaysia's PC2 Score ---

# Combine PCA scores with metadata
df_pca <- data.frame(x = d.pr$x[, 1], y = d.pr$x[, 2]) %>%
  bind_cols(md)

# Calculate centroid (mean PC scores) for each country
df_country_means <- df_pca %>%
  group_by(country, region) %>%
  summarise(
    y.mean = mean(y), 
    x.mean = mean(x), 
    .groups = "drop"
  )

# Get Malaysia's mean PC2 score as the threshold
malaysia_pc2_mean <- df_country_means$y.mean[df_country_means$country == "Malaysia"]

# Define Non-Industrialized countries based on the threshold
country_list_Non_Industrialized <- df_country_means$country[df_country_means$y.mean >= malaysia_pc2_mean]

cat("Threshold (Malaysia PC2):", malaysia_pc2_mean, "\n")
cat("Identified Non-Industrialized Countries:", paste(country_list_Non_Industrialized, collapse=", "), "\n")

# --- 5. Assign Groups to Metadata ---

df_plot <- df_pca %>%
  mutate(
    country_label = case_when(
      country == "Democratic Republic of the Congo" ~ "Dem. Rep. Congo",
      country == "Republic of Congo" ~ "Rep. Congo",
      country == "Central African Republic" ~ "Central African Rep.",
      country == "United States" ~ "USA",
      country == "United Kingdom" ~ "UK",
      TRUE ~ country
    ),
    Analysis_Group = case_when(
    # 1. Sub-Saharan Africa countries
    region == "Sub-Saharan Africa" ~ "Sub-Saharan Africa",
    
    # 2. Non-Industrialized countries (excluding Africa)
    country %in% country_list_Non_Industrialized ~ "Other Non-Industrialized",
    
    # 3. All other countries (Industrialized)
    TRUE ~ "Industrialized"
  ))

# --- 6. Aggregate for Plotting (Centroids and Standard Errors) ---

df_summary <- df_plot %>%
  group_by(country, country_label, Analysis_Group) %>%
  summarise(
    x.mean = mean(x), 
    y.mean = mean(y), 
    n = n(), 
    x.sd = sd(x, na.rm = TRUE), 
    y.sd = sd(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x.se = x.sd / sqrt(n),
    y.se = y.sd / sqrt(n)
  )

group_labels <- df_summary %>%
  group_by(Analysis_Group) %>%
  summarise(
    x.mean = mean(x.mean),
    y.mean = mean(y.mean),
    .groups = "drop"
  )

# --- 7. Plotting ---

p_pca <- ggplot(df_summary, aes(x = x.mean, y = y.mean, fill = Analysis_Group, color = Analysis_Group)) +
  
  stat_ellipse(
    geom = "polygon", 
    alpha = 0.1,        
    linewidth = 0.3,     
    linetype = "dashed",
    show.legend = FALSE
  ) +
  
  geom_errorbar(aes(ymin = y.mean - y.se, ymax = y.mean + y.se), width = 0, linewidth = 0.25) +
  geom_errorbar(aes(xmin = x.mean - x.se, xmax = x.mean + x.se), width = 0, linewidth = 0.25) +
  geom_point(size = 2.5, shape = 21, color = "black", alpha = 0.8, stroke = 0.25) +
  
  ggrepel::geom_text_repel(
    aes(label = country_label), 
    size = 2,            
    color = "black",
    segment.size = 0.15,  
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  
  ggrepel::geom_text_repel(
    data = group_labels,
    aes(label = stringr::str_wrap(Analysis_Group, 15)), 
    size = 2.5,
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.15,      
    box.padding = 1.5,         
    point.padding = 0,
    segment.color = NA,         
    show.legend = FALSE
  ) +
  
  xlab(xlab) + ylab(ylab) +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  # Theme adjustments
  theme_nature(base_size = 6) +
  theme(
    legend.position = "none",
    aspect.ratio = 1
  )

ggsave("results/figures/Figure1b_PCA_Country.pdf", plot = p_pca, 
       width = 120, height = 120, units = "mm", 
       useDingbats = FALSE, bg = "transparent")

ggsave("results/figures/Figure1b_PCA_Country.png", plot = p_pca, 
       width = 120, height = 120, units = "mm", dpi = 300, bg = "white")

# --- 8. Export Updated Metadata for Downstream Analysis ---

# Extract the metadata with the newly assigned 'Analysis_Group'
# and save it as a compressed TSV for the next scripts
md_with_groups <- df_plot %>%
  select(-x, -y) # Remove PCA coordinates to keep metadata clean

write_tsv(md_with_groups, "data/metadata_with_groups_processed.tsv.gz")
cat("Successfully saved updated metadata to data/metadata_with_groups_processed.tsv.gz\n")

