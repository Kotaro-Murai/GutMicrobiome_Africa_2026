# ==============================================================================
# Script: 12_Gastranaerophilales_LMM_Cofactor.R
# Description: LMM Forest Plot of Gastranaerophilales associations (Figure 3c) &
#              Diet correlation Barplots (Figure 3d, Metagenome only).
#              16S rRNA diet correlations are saved as a supplementary figure.
# ==============================================================================

library(tidyverse)
library(data.table)
library(lmerTest)
library(broom.mixed)
library(ggplot2)
library(patchwork)

# Load custom plot styles
source("scripts/00_plot_styles.R")
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# --- Define Colors ---
forest_colors <- c(
  "Physiology" = "#34495E", 
  "Disease" = "#7F8C8D"    
)
tidy_col <- c("Metagenome" = "#18BC9C", "16S rRNA" = "#CCBE93") # Assigned names for consistency

# ==============================================================================
# PART 1: LMM Forest Plot (Figure 3c)
# ==============================================================================
# --- 1. Load Data ---
cat("Loading data for LMM...\n")
d_raw_dt <- data.table::fread("data/metalog_motus3_processed.tsv.gz") 
d_raw <- d_raw_dt %>% data.frame(check.names = FALSE)
rownames(d_raw) <- d_raw[, 1]
d_raw <- d_raw[, -1]

metadata <- read_tsv("data/metadata_with_groups_processed.tsv.gz", show_col_types = FALSE)
md_tax <- data.table::fread("data/mOTUs_3.0.0_GTDB_tax.tsv.gz") %>% data.frame(check.names = FALSE)

# --- 2. GTDB Order level aggregation ---
cat("Aggregating to GTDB Order level...\n")
colnames(d_raw) <- colnames(d_raw) %>% str_remove(".* \\[") %>% str_remove("\\]")
if("-1" %in% colnames(d_raw)) d_raw <- d_raw[, -which(colnames(d_raw) == "-1")]

md2 <- md_tax %>% filter(id %in% colnames(d_raw)) %>% arrange(match(id, colnames(d_raw)))
d_agg <- t(rowsum(t(d_raw), group = md2$order)) %>% data.frame(check.names = FALSE)

target_taxon <- "o__Gastranaerophilales"
if(!target_taxon %in% colnames(d_agg)) stop("Target taxon not found!")

d_target <- d_agg %>%
  select(abundance = all_of(target_taxon)) %>%
  rownames_to_column("rowname") %>%
  mutate(log10_abundance = log10(abundance + 1E-4))

# --- 3. Prepare analysis dataframe ---
cat("Preparing analysis dataframe...\n")

analysis_df <- metadata %>%
  select(rowname = sample_alias, country, region, study, disease, age, sex, bmi, Analysis_Group) %>%
  inner_join(d_target, by = "rowname") %>%
  mutate(
    disease = replace_na(disease, "Control"),
    
    # Standardize continuous variables (Z-score)
    age_z = as.numeric(scale(as.numeric(age))[,1]),
    bmi_z = as.numeric(scale(as.numeric(bmi))[,1]),
    
    sex = factor(sex),
    geo_group = factor(Analysis_Group, levels = c("Industrialized", "Other Non-Industrialized", "Sub-Saharan Africa"))
  )

disease_counts <- analysis_df %>% count(disease)
major_diseases <- disease_counts %>% filter(n >= 100 | disease == "Control") %>% pull(disease)

analysis_df <- analysis_df %>%
  mutate(
    disease_group = if_else(disease %in% major_diseases, disease, "Other"),
    disease_group = factor(disease_group),
    disease_group = relevel(disease_group, ref = "Control") # Reference is Control
  )

# --- 4. Run Separate LMMs ---
cat("Running LMMs...\n")

run_lmm_wrapper <- function(formula_str, data, suffix) {
  tryCatch({
    model <- lmer(as.formula(formula_str), data = data)
    res <- tidy(model, effects = "fixed", conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(model = suffix, .before = 1)
    return(res)
  }, error = function(e) { 
    cat("Error in model", suffix, ":", e$message, "\n")
    return(NULL) 
  })
}

# Model 1: Disease & Geography
res_disease <- run_lmm_wrapper(
  "log10_abundance ~ disease_group + geo_group + (1|study) + (1|country)", 
  data = analysis_df, 
  suffix = "Disease"
)

# Model 2: Physiology & Geography
df_physio <- analysis_df %>% filter(!is.na(bmi_z) & !is.na(age_z) & !is.na(sex))
res_physio <- run_lmm_wrapper(
  "log10_abundance ~ age_z + sex + bmi_z + geo_group + (1|study) + (1|country)", 
  data = df_physio, 
  suffix = "Physio"
)

# --- 5. Format for Forest Plot ---
plot_df <- bind_rows(
  res_disease %>% 
    filter(str_detect(term, "disease")) %>% 
    filter(term != "disease_groupOther") %>%
    mutate(category = "Disease", FDR = p.adjust(p.value, method = "BH")),
  res_physio %>% 
    filter(!str_detect(term, "geo_group|disease")) %>% 
    mutate(category = "Physiology", FDR = p.adjust(p.value, method = "BH"))
) %>%
  mutate(
    term_clean = case_when(
      term == "age_z" ~ "Age (per 1 SD)",
      term == "bmi_z" ~ "BMI (per 1 SD)",
      str_detect(term, "sexmale") ~ "Sex (Male vs Female)",
      str_detect(term, "disease_group") ~ str_replace(term, "disease_group", ""),
      TRUE ~ term
    ),
    ast = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    category = factor(category, levels = c("Physiology", "Disease")),
    plot_color = category
  ) %>%
  arrange(category, estimate) %>%
  mutate(term_clean = factor(term_clean, levels = unique(term_clean)))

# --- 6. Plotting Forest Plot ---
p_forest <- ggplot(plot_df, aes(x = term_clean, y = estimate, color = plot_color)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.25) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, linewidth = 0.25) +
  geom_point(size = 1.5) +
  geom_text(aes(label = ast, 
                y = ifelse(estimate > 0, conf.high + 0.05, conf.low - 0.05)), 
            vjust = ifelse(plot_df$estimate > 0, 0, 1), 
            hjust = 0.5, size = 2.5, color = "black", show.legend = FALSE) +
  
  facet_grid(. ~ category, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_color_manual(values = forest_colors) +
  coord_cartesian(clip = "off") +
  
  labs(title = NULL, y = expression(Log[10] ~ "LMM Estimate"), x = NULL) +
  theme_nature() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text.x = element_text(margin = margin(t = 2, b = 2)),
    strip.background = element_rect(fill = "gray90", color = "gray60", linewidth = 0.25),
    panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.25),
    panel.grid.major.y = element_line(color = "gray80", linetype = "dotted", linewidth = 0.15),
    panel.spacing = unit(2, "mm")
  )

# ==============================================================================
# PART 2: FAOSTAT Food Correlation Barplots (Figure 3d & Supplementary)
# ==============================================================================
cat("Calculating FAOSTAT diet correlations...\n")

## Function to calculate correlation
calc_faostat_cor <- function(d){
  res <- list()
  
  md <- read.delim("data/FAOSTAT.food_g_capita_day.11-14-2025.tsv", row.names = 1, check.names = F) %>% t() %>% data.frame(check.names = F)
  rownames(md) <- rownames(md) %>% str_replace_all("\\.", " ") %>% str_squish()
  rownames(md) <- rownames(md) %>% 
    str_replace("Côte d Ivoire", "Cote d'Ivoire") %>% 
    str_replace("Guinea Bissau", "Guinea-Bissau") %>%
    str_replace("Congo", "Republic of Congo")
  
  keep1 <- rownames(md) %in% d$country
  keep2 <- d$country %in% rownames(md)
  
  md <- md[keep1, ]
  d <- d[keep2, ]
  
  md <- md[order(rownames(md)), ]
  d <- d[order(d$country), ]
  
  ## Spearman correlation
  res_test <- psych::corr.test(d[, 2], md, method = "spearman")
  df <- data.frame(r = t(res_test$r), p = t(res_test$p)) %>% 
    rownames_to_column("diet") %>% 
    filter(!is.na(r)) %>% 
    arrange(-r)
  
  df$fdr <- df$p %>% p.adjust(method = "fdr")
  df <- df %>%
    mutate(
      ast = case_when(
        fdr < 0.001 ~ "***",
        fdr < 0.01 ~ "**",
        fdr < 0.05 ~ "*",
        fdr < 0.1  ~ "+",
        TRUE       ~ ""
      )
    )
  df$n <- nrow(d)
  
  res[[1]] <- df
  res[[2]] <- md %>% rownames_to_column("country") %>% gather("key", "value", -country) %>% left_join(d, by = "country")
  
  return(res)
}

## Run correlation analysis
d_metag <- read.delim("data/Gastranaerophilales_prevalence.metaG.tsv")
cor.metag <- calc_faostat_cor(d_metag)
df_metag <- cor.metag[[1]] %>% mutate(data = "Metagenome")

d_16s <- read.delim("data/Gastranaerophilales_prevalence.16s.tsv")
d_16s <- d_16s %>% filter(n > 20)
cor.16s <- calc_faostat_cor(d_16s)
df_16s <- cor.16s[[1]] %>% mutate(data = "16S rRNA")

# Identify significant diets from EITHER dataset to keep axes consistent if desired, 
# but here we filter independently for their respective plots.
diet.list.metag <- df_metag %>% filter(ast != "" & r > 0) %>% pull(diet)
diet.list.16s <- df_16s %>% filter(ast != "" & r > 0) %>% pull(diet)

# --- Plot Main Figure 3d (Metagenome ONLY) ---
plot_df_metag <- df_metag %>% filter(diet %in% diet.list.metag)

p_diet_main <- ggplot(plot_df_metag, aes(x = reorder(diet, -r), y = r, fill = data)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25, width = 0.7) +
  geom_text(aes(label = ast), size = 2.5, vjust = -0.5) +
  scale_fill_manual(values = tidy_col["Metagenome"]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(y = "Spearman's correlation", x = NULL) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none", # Removed legend since it's only one group
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

# --- Plot Supplementary Figure (16S rRNA ONLY) ---
plot_df_16s <- df_16s %>% filter(diet %in% diet.list.16s)

p_diet_supp <- ggplot(plot_df_16s, aes(x = reorder(diet, -r), y = r, fill = data)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25, width = 0.7) +
  geom_text(aes(label = ast), size = 2.5, vjust = -0.5) +
  scale_fill_manual(values = tidy_col["16S rRNA"]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(y = "Spearman's correlation", x = NULL, title = "16S rRNA Validation") +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    legend.position = "none"
  )

ggsave("results/figures/Supplementary_Figure1_Gastranaerophilales_diet_16s.pdf", plot = p_diet_supp, width = 60, height = 70, units = "mm")


# ==============================================================================
# PART 3: Combine Plots (Figure 3c and 3d)
# ==============================================================================
cat("Combining Figure 3c and Figure 3d...\n")

combined_plot <- p_forest + p_diet_main + 
  plot_layout(widths = c(2, 1))

ggsave("results/figures/Figure3_cd_combined.pdf", 
       plot = combined_plot, 
       width = 180, 
       height = 70, 
       units = "mm")
ggsave("results/figures/Figure3_cd_combined.png", 
       plot = combined_plot, 
       width = 180, 
       height = 70, 
       units = "mm", dpi = 300, bg = "white")

cat("Successfully saved Combined Figure 3cd and Supplementary Figure!\n")

# ==============================================================================
# PART 4: Supplementary Table 7 (LMM Results for Figure 3c)
# ==============================================================================
cat("Generating Supplementary Table 7 (LMM Results)...\n")

supp_table_7 <- plot_df %>%
  select(
    Model_Category = category,
    Variable = term_clean,
    Estimate = estimate,
    Conf_Low = conf.low,
    Conf_High = conf.high,
    Std_Error = std.error,
    T_statistic = statistic,
    P_value = p.value,
    FDR = FDR,
    Significance = ast
  ) %>%
  mutate(
    Estimate = round(Estimate, 4),
    Conf_Low = round(Conf_Low, 4),
    Conf_High = round(Conf_High, 4),
    Std_Error = round(Std_Error, 4)
  ) %>%
  arrange(Model_Category, desc(Estimate))

write_csv(supp_table_7, "results/tables/Supplementary_Table_7_LMM_Gastranaerophilales.csv")

cat("Supplementary Table 7 saved successfully!\n")

# ==============================================================================
# PART 5: Supplementary Table 8 (Diet Correlations: Metagenome & 16S)
# ==============================================================================
cat("Generating Supplementary Table 8 (Diet Correlations)...\n")

supp_table_8 <- full_join(
  df_metag %>% 
    select(Diet = diet, 
           r_Metagenome = r, 
           P_value_Metagenome = p, 
           FDR_Metagenome = fdr, 
           N_Metagenome = n, 
           Significance_Metagenome = ast),
  
  df_16s %>% 
    select(Diet = diet, 
           r_16S = r, 
           P_value_16S = p, 
           FDR_16S = fdr, 
           N_16S = n, 
           Significance_16S = ast),
  
  by = "Diet"
) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 4))
  ) %>%
  arrange(desc(r_Metagenome))

write_csv(supp_table_8, "results/tables/Supplementary_Table_8_Diet_Correlations.csv")

cat("Supplementary Table 8 saved successfully!\n")