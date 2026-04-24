# ==============================================================================
# Script: Supple_lifestyle_location_prevalence.R
# Description: Combined Barplot of Gastranaerophilales prevalence by Lifestyle 
#              (Left) and Location faceted by Lifestyle (Right). (Supplementary Figure 1)
# ==============================================================================

library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork) 
library(countrycode)

# Load custom plot styles
source("scripts/00_plot_styles.R")
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- 1. Load Data ---
cat("Loading data...\n")
md <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)
d_raw_dt <- data.table::fread("data/metalog_motus3_processed.tsv.gz") 
d_raw <- d_raw_dt %>% data.frame(check.names = FALSE)
rownames(d_raw) <- d_raw[, 1]
d_raw <- d_raw[, -1]

md_tax <- data.table::fread("data/mOTUs_3.0.0_GTDB_tax.tsv.gz") %>% data.frame(check.names = FALSE)

# --- 2. Extract Target Taxon ---
cat("Extracting and aggregating target taxon...\n")
colnames(d_raw) <- if_else(
  str_detect(colnames(d_raw), "\\["),
  str_extract(colnames(d_raw), "(?<=\\[)[^\\[\\]]+(?=\\]$)"),
  colnames(d_raw)
)
if("-1" %in% colnames(d_raw)) d_raw <- d_raw[, -which(colnames(d_raw) == "-1")]

md2 <- md_tax %>% filter(id %in% colnames(d_raw)) %>% arrange(match(id, colnames(d_raw)))
d_agg <- t(rowsum(t(d_raw), group = md2$order)) %>% data.frame(check.names = FALSE)

target_taxon <- "o__Gastranaerophilales"
d_target <- d_agg %>%
  select(Raw_Abundance = all_of(target_taxon)) %>%
  rownames_to_column("sample_alias")

# --- 3. Merge and Process Metadata ---
cat("Merging with metadata and creating formatted labels...\n")
df_analysis <- d_target %>%
  left_join(md, by = "sample_alias") %>% 
  filter(!is.na(lifestyle)) %>% 
  mutate(
    country_iso3 = countrycode(country, origin = "country.name", destination = "iso3c"),
    study_short = str_extract(study, "^.*?_\\d{4}"),
    plot_label = paste0(location, "_", country_iso3, " (", study_short, ")"),
    Presence = if_else(Raw_Abundance > 0.0001, 1, 0)
  ) %>%
  mutate(lifestyle = factor(lifestyle, levels = c(
    "Hunter-gatherer", 
    "Rural_agrarian / traditional", 
    "Peri-urban / semi-urban", 
    "Urban"
  )))
label_order <- df_analysis %>%
  group_by(plot_label) %>%
  summarize(prevalence = mean(Presence) * 100) %>%
  arrange(desc(prevalence)) %>%
  pull(plot_label)

df_analysis <- df_analysis %>%
  mutate(plot_label = factor(plot_label, levels = label_order))

# ==============================================================================
# 4. Calculate Prevalence
# ==============================================================================
cat("Calculating prevalence...\n")
prev_lifestyle <- df_analysis %>%
  group_by(lifestyle) %>%
  summarise(
    n_total = n(),
    n_present = sum(Presence),
    prevalence = (n_present / n_total) * 100,
    .groups = "drop"
  )
prev_location <- df_analysis %>%
  group_by(lifestyle, plot_label) %>%
  summarise(
    n_total = n(),
    n_present = sum(Presence),
    prevalence = (n_present / n_total) * 100,
    .groups = "drop"
  )

# ==============================================================================
# 5. Plotting: Left (Lifestyle Overall) with Statistical Significance
# ==============================================================================
cat("Calculating pairwise statistics for prevalence (Fisher's exact test)...\n")

lifestyle_levels <- levels(df_analysis$lifestyle)
comparisons <- combn(lifestyle_levels, 2, simplify = FALSE)

stat_pvals <- map_dfr(comparisons, function(comp) {
  sub_df <- df_analysis %>% filter(lifestyle %in% comp)
  tbl <- table(sub_df$lifestyle, sub_df$Presence)
  tbl <- tbl[rowSums(tbl) > 0, ] 
  res <- fisher.test(tbl)
  
  tibble(
    group1 = comp[1],
    group2 = comp[2],
    p = res$p.value
  )
}) %>%
  mutate(
    p.adj = p.adjust(p, method = "bonferroni"),
    p.adj.signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001  ~ "***",
      p.adj <= 0.01   ~ "**",
      p.adj <= 0.05   ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  filter(p.adj.signif != "ns")
if (nrow(stat_pvals) > 0) {
  stat_pvals <- stat_pvals %>%
    mutate(y.position = 105 + (row_number() - 1) * 12) 
  max_y_left <- max(stat_pvals$y.position) + 5
} else {
  max_y_left <- 100
}

cat("Generating Left Plot (Lifestyle Prevalence)...\n")

p_left <- ggplot(prev_lifestyle, aes(x = lifestyle, y = prevalence, fill = lifestyle)) +
  geom_col(alpha = 0.8, color = "gray20", linewidth = 0.25) +
  {if (nrow(stat_pvals) > 0) stat_pvalue_manual(
    stat_pvals,
    label = "p.adj.signif",
    y.position = "y.position",
    tip.length = 0.02,
    size = 2,
    vjust = 0.5
  )} +
  scale_y_continuous(
    limits = c(0, max_y_left), 
    breaks = seq(0, 100, by = 25), 
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = lifestyle_colors) + 
  labs(
    y = "Prevalence (%)",
    x = NULL
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# ==============================================================================
# 6. Plotting: Right (Location Prevalence)
# ==============================================================================
cat("Generating Right Plot (Location Prevalence)...\n")

p_right <- ggplot(prev_location, aes(x = plot_label, y = prevalence, fill = lifestyle)) +
  geom_col(alpha = 0.8, color = "gray20", linewidth = 0.25) +
  facet_grid(~ lifestyle, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    limits = c(0, max_y_left), 
    breaks = seq(0, 100, by = 25), 
    expand = expansion(mult = c(0, 0.05))
  ) +
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
# 7. Combine and Save
# ==============================================================================
cat("Combining and saving Supplementary Figure 1...\n")

p_combined <- p_left + p_right + 
  plot_layout(widths = c(1, 5))

outfile_plot <- "results/figures/Supplementary_Figure1_Prevalence"

ggsave(paste0(outfile_plot, ".pdf"), p_combined, width = 180, height = 100, units = "mm")
ggsave(paste0(outfile_plot, ".png"), p_combined, width = 180, height = 100, units = "mm", dpi = 300, bg = "white")
