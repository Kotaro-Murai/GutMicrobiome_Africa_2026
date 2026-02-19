# ==============================================================================
# Script: 02_pca_clustering.R
# Author: Kotaro Murai
# Description: 
#   Perform Principal Component Analysis (PCA) on species-level relative abundances.
#   Define three analysis groups (Industrialized, Non-Africa Non-Industrialized, 
#   and Africa) based on the PC2 score threshold (Malaysia).
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
  mutate(Analysis_Group = case_when(
    # 1. Sub-Saharan Africa countries
    region == "Sub-Saharan Africa" ~ "Africa",
    
    # 2. Non-Industrialized countries (excluding Africa)
    country %in% country_list_Non_Industrialized ~ "Non-Africa Non-Industrialized",
    
    # 3. All other countries (Industrialized)
    TRUE ~ "Industrialized"
  ))

# Set factor levels for consistent legend ordering
df_plot$Analysis_Group <- factor(
  df_plot$Analysis_Group, 
  levels = c("Industrialized", "Non-Africa Non-Industrialized", "Africa")
)

# --- 6. Aggregate for Plotting (Centroids and Standard Errors) ---

df_summary <- df_plot %>%
  group_by(country, Analysis_Group) %>%
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

# --- 7. Plotting ---

p_pca <- ggplot(df_summary, aes(x = x.mean, y = y.mean, fill = Analysis_Group, color = Analysis_Group)) +
  theme_thesis(base_size = 14) + 
  xlab(xlab) + ylab(ylab) +
  
  # Error bars
  geom_errorbar(aes(ymin = y.mean - y.se, ymax = y.mean + y.se), width = 0, size = 0.5) +
  geom_errorbar(aes(xmin = x.mean - x.se, xmax = x.mean + x.se), width = 0, size = 0.5) +
  
  # Country centroids
  geom_point(aes(size = n), shape = 21, color = "black", alpha = 0.8) +
  
  # Aesthetics
  scale_fill_manual(values = thesis_colors) +
  scale_color_manual(values = thesis_colors) +
  scale_size_continuous(trans = "log10", name = "Sample size", range = c(1, 4)) +
  
  # Country labels
  ggrepel::geom_text_repel(
    aes(label = country), 
    size = 4, 
    color = "black",
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  
  # Theme adjustments
  theme(
    legend.position = "bottom"
  )

# Display and save plot
print(p_pca)
ggsave("results/figures/Figure1b_PCA_Clustering.pdf", plot = p_pca, width = 10, height = 7, useDingbats = FALSE)
ggsave("results/figures/Figure1b_PCA_Clustering.png", plot = p_pca, width = 10, height = 7, dpi = 300)

# --- 8. Export Updated Metadata for Downstream Analysis ---

# Extract the metadata with the newly assigned 'Analysis_Group'
# and save it as a compressed TSV for the next scripts
md_with_groups <- df_plot %>%
  select(-x, -y) # Remove PCA coordinates to keep metadata clean

write_tsv(md_with_groups, "data/metadata_with_groups_processed.tsv.gz")
cat("Successfully saved updated metadata to data/metadata_with_groups_processed.tsv.gz\n")