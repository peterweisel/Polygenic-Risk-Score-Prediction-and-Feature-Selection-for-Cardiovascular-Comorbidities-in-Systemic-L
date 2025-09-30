#!/bin/bash -l

# Author: Peter Joseph Weisel
# Affiliation: Forskargruppen reumatologi @ Institutionen för medicinska vetenskaper
# Description: This script uses PLINK to perform allelic scoring for individuals using a set of SNPs associated with Coronary artery disease
# Date: 2025-09-29

# Export executable to PATH
cp /Users/peterweisel/Desktop/CVD/src/plink_mac_20250819/plink ~/bin/
export PATH="$HOME/bin:$PATH"

# Directories
PLINK=/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_plink_dataset
DATA_PLINK=/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_plink_dataset
DATA_SCORE=/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_edited_SCOREs

# PGS trait
PHEN="CAD"

# List of rsIDs collected from PGS Catalog
declare -a lst=("PGS004199" "PGS004307")
for SAMPLE in ${lst[@]};
do
    OUT=/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_plink_dataset/assigned_scores/${PHEN}/${SAMPLE}
    mkdir -p ${OUT}
    cd ${OUT}
    
    # Allelic scoring
    plink --bfile ${DATA_PLINK}/sle_snv_snvs_set1 --score ${DATA_SCORE}/_ed_${SAMPLE}.txt --out ${PHEN}_${SAMPLE}_indiv_score
done
