library(tidyverse)

source("scripts/00_plot_styles.R")

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

tidy_col <- c("#18BC9C","#CCBE93","#a6cee3","#1f78b4","#b2df8a","#fb9a99","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")

## association with FAOSTAT data
calc_faostat_cor <- function(d){
  res <- list()
  
  md <- read.delim("data/FAOSTAT.food_g_capita_day.11-14-2025.tsv", row.names = 1, check.names = F) %>% t() %>% data.frame(check.names = F)
  rownames(md) <- rownames(md) %>% str_replace_all("\\.", " ") %>%
    str_squish()
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
  res <- psych::corr.test(d[, 2], md, method = "spearman")
  df <- data.frame(r = t(res$r), p = t(res$p)) %>% 
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

## Metagenome
d <- read.delim("data/Gastranaerophilales_prevalence.metaG.tsv")
cor.metag <- calc_faostat_cor(d)

## 16S rRNA
d <- read.delim("data/Gastranaerophilales_prevalence.16s.tsv")
d <- d %>% filter(n > 20) ## exclude countries less than 20
cor.16s <- calc_faostat_cor(d)

cor.metag[[1]]$data <- "Metagenome"
cor.metag[[2]]$data <- "Metagenome"
cor.16s[[1]]$data <- "16S rRNA"
cor.16s[[2]]$data <- "16S rRNA"

df <- rbind(cor.metag[[1]], cor.16s[[1]])
diet.list <- df %>% filter(ast != "" & r > 0) %>% pull(diet)
df <- df %>% filter(diet %in% diet.list)

df$data <- df$data %>% factor(levels = c("Metagenome", "16S rRNA"))

p <- ggplot(df, aes(x = reorder(diet, -r, sum), y = r, fill = data)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.25) +
  geom_text(aes(label = ast), position = position_dodge(width = 0.9), size = 2.5, vjust = -0.5) +
  scale_fill_manual(values = tidy_col) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    y = "Spearman's correlation", 
    x = NULL
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = c(0.95, 0.95),
    legend.justification = c(0.95, 0.95),
    legend.title = element_blank(),
    legend.margin = margin(b = -5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

ggsave("results/figures/Figure3d_Gastranaerophilales_diet_correlation_barplots.pdf", plot = p, width = 60, height = 70, units = "mm")
