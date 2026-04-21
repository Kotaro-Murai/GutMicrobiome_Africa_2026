# Gut Microbiome Africa 2026

This repository contains the R scripts used for the statistical analysis and visualization in the "Gut Microbiome in Africa" study. 

## Overview
This repository provides:
- R scripts implementing all statistical analyses described in the study.
- Figure- and table-generating scripts corresponding to the manuscript.
- A transparent and reproducible analysis workflow.

## Repository Structure
```text
.
├── data/           
├── results/        # Generated outputs
│   ├── figures/    # Main and supplementary figures
│   └── tables/     # Supplementary tables
├── scripts/
│   ├── 00_disease_data_process.R
│   ├── 00_prepare_16s_data.R
│   ├── 00_plot_styles.R
│   ├── 01_map_metagenome_data.R
│   ├── 02_pca_country.R
│   ├── 03_alpha_diversity.R
│   ├── 04_taxa_composition.R
│   ├── 05_lmm_species_volcano.R
│   ├── 06_Lollipop_Species.R
│   ├── 07_phylogenetic_tree_barplot.R
│   ├── 08_lmm_higher_ranks.R
│   ├── 09_Metalog_Map_Gastranaerophilales.R
│   ├── 10_16s_Map_Gastranaerophilales.R
│   ├── 11_lifestyle_location.R
│   ├── 12_Gastranaerophilales_LMM_Cofactor.R
│   ├── 13_Food_Correlation.R
│   ├── 14_pcoa_cyanobacteria_uhgg.R
│   ├── 15_functional_enrichment_gastranaerophilales.R
│   ├── 16_final_figures.R    
│   ├── Supple_target_genera_boxplots.R
│   ├── Supple_lifestyle_location_prevalence.R
│   ├── Supple_plot_flagellar_heatmap.R      
│   └── Supple_phylogenetic_distribution.R          
└── README.md
