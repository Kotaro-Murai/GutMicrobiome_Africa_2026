# ==============================================================================
# Script: 11_lifestyle_location.R
# Description: Combined Boxplot of Gastranaerophilales abundance by Lifestyle 
#              (Left) and Location faceted by Lifestyle (Right). (Figure 3b)
# ==============================================================================

library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork) 
library(countrycode)

# Load custom plot styles
source("scripts/00_plot_styles.R")
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---
cat("Loading data...\n")

md <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)

if(!file.exists("data/metalog_motus3_processed.tsv.gz")) stop("Error: Abundance file not found.")
d_raw_dt <- data.table::fread("data/metalog_motus3_processed.tsv.gz") 
d_raw <- d_raw_dt %>% data.frame(check.names = FALSE)
rownames(d_raw) <- d_raw[, 1]
d_raw <- d_raw[, -1]

if(!file.exists("data/mOTUs_3.0.0_GTDB_tax.tsv.gz")) stop("Error: Taxonomy file not found.")
md_tax <- data.table::fread("data/mOTUs_3.0.0_GTDB_tax.tsv.gz") %>% data.frame(check.names = FALSE)

# --- 2. Extract Target Taxon ---
cat("Extracting and aggregating target taxon...\n")

colnames(d_raw) <- colnames(d_raw) %>% str_remove(".* \\[") %>% str_remove("\\]")
if("-1" %in% colnames(d_raw)) d_raw <- d_raw[, -which(colnames(d_raw) == "-1")]

md2 <- md_tax %>% filter(id %in% colnames(d_raw)) %>% arrange(match(id, colnames(d_raw)))
d_agg <- t(rowsum(t(d_raw), group = md2$order)) %>% data.frame(check.names = FALSE)

target_taxon <- "o__Gastranaerophilales"
if(!target_taxon %in% colnames(d_agg)) stop("Target taxon not found!")

d_target <- d_agg %>%
  select(Raw_Abundance = all_of(target_taxon)) %>%
  rownames_to_column("sample_alias") %>%
  mutate(log10_Abundance = log10(Raw_Abundance + 1e-4))

# --- 3. Merge and Process Metadata ---
cat("Merging with metadata and creating formatted labels...\n")

df_analysis <- d_target %>%
  left_join(md, by = "sample_alias") %>% 
  filter(!is.na(lifestyle)) %>% 
  mutate(
    country_iso3 = countrycode(country, origin = "country.name", destination = "iso3c"),
    study_short = str_extract(study, "^.*?_\\d{4}"),
    plot_label = paste0(location, "_", country_iso3, " (", study_short, ")")
  ) %>%
  mutate(lifestyle = factor(lifestyle, levels = c(
    "Hunter-gatherer", 
    "Rural_agrarian / traditional", 
    "Peri-urban / semi-urban", 
    "Urban"
  ))) %>%
  mutate(plot_label = fct_reorder(plot_label, log10_Abundance, .fun = median, .desc = TRUE)) %>%
  select(sample_alias, study, location, country, country_iso3, plot_label, lifestyle, Raw_Abundance, log10_Abundance)

y_min <- min(df_analysis$log10_Abundance, na.rm = TRUE) - 0.2
y_max_data <- max(df_analysis$log10_Abundance, na.rm = TRUE) 


# ==============================================================================
# 4. Left: Lifestyle Overall & Stats Calculation
# ==============================================================================
cat("Calculating pairwise statistics for the plot...\n")

compare_res <- compare_means(
  log10_Abundance ~ lifestyle, 
  data = df_analysis, 
  method = "wilcox.test", 
  p.adjust.method = "bonferroni"
) %>%
  mutate(
    p.adj.signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001  ~ "***",
      p.adj <= 0.01   ~ "**",
      p.adj <= 0.05   ~ "*",
      TRUE            ~ "ns"
    )
  )

y_step <- 0.4 

stat_pvals <- compare_res %>%
  filter(p.adj.signif != "ns") %>%
  mutate(y.position = y_max_data + row_number() * y_step)

y_max_plot <- ifelse(nrow(stat_pvals) > 0, max(stat_pvals$y.position) + 0.3, y_max_data + 0.5)
y_limits <- c(y_min, y_max_plot) 

cat("Generating Left Plot (Lifestyle Overall)...\n")

p_left <- ggplot(df_analysis, aes(x = lifestyle, y = log10_Abundance, fill = lifestyle)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "gray20", linewidth = 0.25) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.5, color = "gray40", stroke = 0) +
  
  stat_pvalue_manual(
    stat_pvals,
    label = "p.adj.signif",   
    y.position = "y.position", 
    tip.length = 0.02,
    size = 2,                
    vjust = 0.5
  ) +
  
  scale_y_continuous(limits = y_limits) +
  scale_fill_manual(values = lifestyle_colors) + 
  
  labs(
    y = expression(Log[10] ~ "Relative Abundance"),
    x = NULL
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# ==============================================================================
# 5. Right: Location by Lifestyle
# ==============================================================================
cat("Generating Right Plot (Location by Lifestyle)...\n")

p_right <- ggplot(df_analysis, aes(x = plot_label, y = log10_Abundance, fill = lifestyle)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "gray20", linewidth = 0.25) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.5, color = "gray40", stroke = 0) +
  
  facet_grid(~ lifestyle, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = y_limits) +
  scale_fill_manual(values = lifestyle_colors) + 
  
  labs(
    y = NULL, 
    x = NULL
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    
    strip.text = element_text(margin = margin(t = 2, b = 2)),
    strip.background = element_rect(fill = "gray90", color = "gray60", linewidth = 0.25),
    panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.25),
    panel.spacing = unit(1, "mm") 
  )

# ==============================================================================
# 6. Combine and Save (Patchwork)
# ==============================================================================
cat("Combining and saving Figure 3b...\n")

p_combined <- p_left + p_right + 
  plot_layout(widths = c(1, 5))

outfile_plot <- "results/figures/Figure3b_Lifestyle_Location"

ggsave(paste0(outfile_plot, ".pdf"), p_combined, width = 180, height = 70, units = "mm")
ggsave(paste0(outfile_plot, ".png"), p_combined, width = 180, height = 70, units = "mm", dpi = 300, bg = "white")

cat("Successfully saved Figure 3b\n")

# ==============================================================================
# 8. Generate Supplementary Table 5 (African Samples Lifestyle Categorization)
# ==============================================================================
cat("Generating Supplementary Table 5...\n")

supp_table_5 <- df_analysis %>%
  select(
    Sample_ID = sample_alias,
    Study = study,
    Country = country,
    Location = location,
    Lifestyle_Category = lifestyle,
    Gastranaerophilales_Relative_Abundance = Raw_Abundance,
    Log10_Abundance = log10_Abundance
  ) %>%
  arrange(Lifestyle_Category, Country, Study, Location)

write_csv(supp_table_5, "results/tables/Supplementary_Table_5_African_Samples_Lifestyle.csv")

cat("Supplementary Table 5 saved to results/tables/ directory.\n")