#!/bin/bash -l

# Author: Peter Joseph Weisel
# Affiliation: Forskargruppen reumatologi @ Institutionen för medicinska vetenskaper
# Description: This script uses PLINK to perform allelic scoring for individuals using a set of SNPs associated with Coronary artery disease
# Date: 2025-09-29

# Export executable to PATH
cp /Users/peterweisel/Desktop/CVD/src/plink_mac_20250819/plink ~/bin/
export PATH="$HOME/bin:$PATH"

# Directories
PLINK_PREFIX="/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_plink_dataset/plink_files/CORRECTED_dataset/_corrected_set1"
rsIDs="/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_edited_SCOREs/rsid_list.txt"
OUT="/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_plink_dataset/extracted"
mkdir -p ${OUT}
cd ${OUT}

plink \
  --bfile ${PLINK_PREFIX} \
  --extract ${rsIDs} \
  --recode A \
  --out all_rsids_genotypes
