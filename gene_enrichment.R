#!/usr/bin/env Rscript
# Author: Peter Joseph Weisel
# Affiliation: Forskargruppen reumatologi @ Institutionen för medicinska vetenskaper
# Description: Gene enrichment
# Date: 2025-10-21

## load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)

# load ICVL FUMA GWAS output file
ICVL <- read.delim("ICVL_FUMA/GS.txt", header = TRUE, fill = TRUE, quote = "", comment.char = "")

ICVLP_grouped <- ICVL %>%
  group_by(Category) %>%
  summarise(
    n_sets = n(),
    min_p = min(p),
    mean_p = mean(p),
    top_geneset = GeneSet[which.min(p)]
  ) %>%
  arrange(min_p)

# subset GWAS stats
ICVLP_GWAS <- ICVL %>%
  filter(Category == "GWAScatalog") %>%
  slice_max(order_by = N_genes, n = 10) %>%
  arrange(desc(N_genes))

# subset GO BP stats
ICVLP_GO_bp <- ICVL %>%
  filter(Category == "GO_bp") %>%
  slice_max(order_by = N_genes, n = 10) %>%
  arrange(desc(N_genes))

# subset Reactome stats
ICVLP_Reactome <- ICVL %>%
  filter(Category == "Reactome") %>%
  slice_max(order_by = N_genes, n = 10) %>%
  arrange(desc(N_genes))

# plot GWAS enrichment
GWAS <- ggplot(ICVLP_GWAS, aes(x = reorder(GeneSet, N_genes), y = N_genes)) +
  geom_col(fill = "#EE6549") +
  coord_flip() +
  labs(
    x = "GWAS Trait",
    y = "Number of Genes",
    title = "Top ICVL GWAS Catalog Gene Sets"
  ) +
  theme_minimal(base_size = 14)
GWAS

# plot GO BP enrichment
GO_bp <- ggplot(ICVLP_GO_bp, aes(x = reorder(GeneSet, N_genes), y = N_genes)) +
  geom_col(fill = "#8A5176") +
  coord_flip() +
  labs(
    x = "GO Biological Process",
    y = "Number of Genes",
    title = "Top ICVL GO Catalog Gene Sets"
  ) +
  theme_minimal(base_size = 14)
GO_bp

# plot reactome enrichment
reactome <- ggplot(ICVLP_Reactome, aes(x = reorder(GeneSet, N_genes), y = N_genes)) +
  geom_col(fill = "#426C6B") +
  coord_flip() +
  labs(
    x = "Reactome Trait",
    y = "Number of Genes",
    title = "Top ICVL Reactome Catalog Gene Sets"
  ) +
  theme_minimal(base_size = 14)
reactome
