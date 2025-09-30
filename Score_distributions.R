# !/usr/bin/Rscript
# Author: Peter Joseph Weisel
# Affiliation: Forskargruppen reumatologi @ Institutionen för medicinska vetenskaper
# Description: Extract allelic scoring for PGS traits, analyze distributions, perform non-parametric tests on case and controls
# Date: 2025-09-29

######### Load packages

library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)

######## Functions

# 1. pheno_plot
pheno_plot <- function(df1, df2, labels = c("Dataset1", "Dataset2"), title = "Title") {
  
  # Input: df1 and df2 (dataframe of PLINK allelic scoring outputs for two PGS Catalog traits),
  #        labels (vector containing 2 PGS ID strings), title (string for figure title).
  # Output: 2 density plots showing case vs control distributions for each PGS ID.
  
  df1_plot <- ggplot(df1, aes(x = SCORE, color = PHENO, fill = PHENO)) +
    geom_density(alpha = 0.4) +
    theme_minimal() +
    labs(x = "Score", y = "Density", color = "Phenotype", fill = "Phenotype")
  
  df2_plot <- ggplot(df2, aes(x = SCORE, color = PHENO, fill = PHENO)) +
    geom_density(alpha = 0.4) +
    theme_minimal() +
    labs(x = "Score", y = "Density", color = "Phenotype", fill = "Phenotype")
  
  combined_plot <- plot_grid(df1_plot, df2_plot, labels = labels)
  
  title_plot <- ggdraw() + 
    draw_label(title, fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  plot_grid(title_plot, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
}

# 2. clean_pheno
clean_pheno <- function(df) {
  
  # Input: df (dataframe of PLINK allelic scoring outputs for PGS Catalog trait).
  # Output: New dataframe categorising individuals as case or controls.
  
  df %>%
    mutate(
      PHENO = factor(PHENO,
                     levels = c(1, 2),
                     labels = c("control", "case"))
    ) %>%
    filter(!is.na(PHENO))
}

######## Import data

# Systemic lupus erythematosus (SLE)
PGS000196 <- read.table("SLE/PGS000196/SLE_PGS000196_indiv_score.profile", header = TRUE)
PGS000803 <- read.table("SLE/PGS000803/SLE_PGS000803_indiv_score.profile", header = TRUE)
PGS003960 <- read.table("SLE/PGS003960/SLE_PGS003960_indiv_score.profile", header = TRUE)
# Atrial fibrillation (AF)
PGS004186 <- read.table("AF/PGS004186/AF_PGS004186_indiv_score.profile", header = TRUE)
PGS004292 <- read.table("AF/PGS004292/AF_PGS004292_indiv_score.profile", header = TRUE)
# Myocardial infarction (MI)
PGS001315 <- read.table("MI/PGS001315/MI_PGS001315_indiv_score.profile", header = TRUE)
PGS000710 <- read.table("MI/PGS000710/MI_PGS000710_indiv_score.profile", header = TRUE)
# Stroke (STROKE)
PGS004000 <- read.table("STROKE/PGS004000/STROKE_PGS004000_indiv_score.profile", header = TRUE)
PGS004124 <- read.table("STROKE/PGS004124/STROKE_PGS004124_indiv_score.profile", header = TRUE)
# Coronary artery disease (CAD)
PGS004199 <- read.table("CAD/PGS004199/CAD_PGS004199_indiv_score.profile", header = TRUE)
PGS004307 <- read.table("CAD/PGS004307/CAD_PGS004307_indiv_score.profile", header = TRUE)

######## Clean data

# SLE
PGS000196 <- clean_pheno(PGS000196)
PGS000803 <- clean_pheno(PGS000803)
PGS003960 <- clean_pheno(PGS003960)
# AF
PGS004186 <- clean_pheno(PGS004186)
PGS004292 <- clean_pheno(PGS004292)
# MI
PGS001315 <- clean_pheno(PGS001315)
PGS000710 <- clean_pheno(PGS000710)
# STROKE
PGS004000 <- clean_pheno(PGS004000)
PGS004124 <- clean_pheno(PGS004124)
# CAD
PGS004199 <- clean_pheno(PGS004199)
PGS004307 <- clean_pheno(PGS004307)


######## SLE plot

# NOTE: pheno_plot will not be used for SLE PGS IDs since there are 3 scoring datasets.
PGS000196_plot <- ggplot(PGS000196, aes(x = SCORE, color = PHENO, fill = PHENO)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(x = "Score", y = "Density", color = "Phenotype", fill = "Phenotype")

PGS000803_plot <- ggplot(PGS000803, aes(x = SCORE, color = PHENO, fill = PHENO)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(x = "Score", y = "Density", color = "Phenotype", fill = "Phenotype")

PGS003960_plot <- ggplot(PGS003960, aes(x = SCORE, color = PHENO, fill = PHENO)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(x = "Score", y = "Density", color = "Phenotype", fill = "Phenotype")

# Combine into one plot
SLE <- plot_grid(PGS000196_plot, PGS000803_plot, PGS003960_plot, labels = c("PGS000196", "PGS000803", "PGS003960"))

# SLE plot title
SLE_title <- ggdraw() + 
  draw_label("Distribution of SLE Scores", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(SLE_title, SLE, ncol = 1, rel_heights = c(0.1, 1))

# SLE non-parametric test between case and control
PGS000196_test <- wilcox.test(SCORE ~ PHENO, data = PGS000196)
PGS000803_test <- wilcox.test(SCORE ~ PHENO, data = PGS000803)
PGS003960_test <- wilcox.test(SCORE ~ PHENO, data = PGS003960)

######## Arterial Fibrillation plot

AF <- pheno_plot(PGS004186, PGS004292,
           labels = c("PGS004186", "PGS004292"),
           title = "Arterial Fibrillation Polygenic Score Distributions")

# AF non-parametric test between case and control
PGS004186_test <- wilcox.test(SCORE ~ PHENO, data = PGS004186)
PGS004292_test <- wilcox.test(SCORE ~ PHENO, data = PGS004292)

######## MI

MI <- pheno_plot(PGS001315, PGS000710,
                 labels = c("PGS001315", "PGS000710"),
                 title = "Myocardial Infraction Polygenic Score Distributions")

# SLE non-parametric test between case and control
PGS001315_test <- wilcox.test(SCORE ~ PHENO, data = PGS001315)
PGS000710_test <- wilcox.test(SCORE ~ PHENO, data = PGS000710)

######## STROKE

STROKE <- pheno_plot(PGS004000, PGS004124,
                 labels = c("PGS004000", "PGS004124"),
                 title = "Stroke Polygenic Score Distributions")

# SLE non-parametric test between case and control
PGS004000_test <- wilcox.test(SCORE ~ PHENO, data = PGS004000)
PGS004124_test <- wilcox.test(SCORE ~ PHENO, data = PGS004124)

######## CAD

CAD <- pheno_plot(PGS004199, PGS004307,
                     labels = c("PGS004199", "PGS004307"),
                     title = "CAD Score Distributions")

# SLE non-parametric test between case and control
PGS004199_test <- wilcox.test(SCORE ~ PHENO, data = PGS004199)
PGS004307_test <- wilcox.test(SCORE ~ PHENO, data = PGS004307)
