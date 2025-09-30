# !/usr/bin/Rscript
# Author: Peter Joseph Weisel
# Affiliation: Forskargruppen reumatologi @ Institutionen för medicinska vetenskaper
# Description: Randomly extract Individual IDs (IIDs) from FAM file
# Date: 2025-09-29

######### Load packages
library(dplyr)

######### Load FAM file
fam_file <- read.table("sle_snv_snvs_set1.fam", header = FALSE)

######### Sample IIDs
# Extract IIDs as vector
indiv_id_list <- unique(fam_file$V2)
# Number of IIDs per group
major_n <- round(length(indiv_id_list) * 0.80)
minor_n <- length(indiv_id_list) - major_n
# Sample without replacement for majority and use remaining IIDs for minority
indivID_80 <- sample(indiv_id_list, major_n, replace = FALSE)
indivID_20 <- setdiff(indiv_id_list, indivID_80)
# Confirm all IIDs were used
length(indivID_80) + length(indivID_20) == length(indiv_id_list)
# Confirm that no overlaps exist between vectors
length(intersect(indivID_80, indivID_20))

######### Export IIDs to desktop
writeLines(indivID_80, "/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_plink_dataset/subset_individuals/indivID_80.txt")
writeLines(indivID_20, "/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_plink_dataset/subset_individuals/indivID_20.txt")
