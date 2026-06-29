# =========================================================
# Part 1. Data preprocessing and LD reference construction
# =========================================================
#
# Description:
# This script prepares ancestry-specific LD reference panels
# from the 1000 Genomes Project Phase 3 dataset for
# downstream polygenic risk score analyses using LDpred2.
#
# Requirements:
# - PLINK v2.0
# - R >= 4.3
# - R packages:
#     data.table
#     Matrix
#     pROC
#     tidyverse
#     dplyr
#     docopt
#     rio
#     magrittr
#     bigsnpr
# - wget
#
# Input:
# - 1000 Genomes Phase 3 VCF files
# - Population panel file
# - HapMap3 SNP list
# - Genetic map files
#
# Output:
# - Ancestry-specific PLINK genotype files
# - HapMap3-filtered genotype files
# - LDpred2 LD reference cache
#
# =========================================================
PROJECT_DIR=$(pwd)

# =========================================================
# 1. Set project directories
# =========================================================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -dir)    PROJECT_DIR="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

SCRIPT_DIR=/opt/pipeline/scripts
GWAS_DIR=${PROJECT_DIR}/GWAS
HM3_DIR=${DATA_DIR}/HapMap3
hm3=${HM3_DIR}/hapmap3.snp
GENO_DIR=${DATA_DIR}/Genotype
DATA_DIR=${PROJECT_DIR}/data
Layer1_DIR=${PROJECT_DIR}/Layer1
Layer2_DIR=${PROJECT_DIR}/Layer2
Weight_DIR=${PROJECT_DIR}/SNP_Weight
VCF_DIR=${DATA_DIR}/1KG_Phase3_VCF
GENETIC_MAP_DIR=${DATA_DIR}/omni_genetic_map
LD_REF_DIR=${DATA_DIR}/LD_ref
PHENO_DIR=${DATA_DIR}/Phenotype
LDpred2_DIR=${DATA_DIR}/LDpred2
PRS_DIR=${DATA_DIR}/PRS


# =========================================================
# Generate LDpred2 LD reference cache
# =========================================================
#
# Requirements:
# - Ldpred2LDRefCache.R
# - R >= 4.3
# - bigsnpr
#
# IMPORTANT:
#
# The directory "omni_genetic_map"
# must be located in the working directory.
#
# =========================================================
cd ${DATA_DIR}

for POP in AFR EUR EAS SAS AMR
do
    mkdir -p ${LD_REF_DIR}/${POP}
    Rscript ${SCRIPT_DIR}/Ldpred2LDRefCache.R \
        -p ${LD_REF_DIR}/${POP}/${POP}_chr#.hm3.
done


# =========================================================
# 6. Phenotype and covariates preparation
# =========================================================
#
# Place phenotype and covariates files in:
#
#   ${PHENO_DIR}
#
# Suggested files:
#
# - phenotype.txt
# - covariates.txt
# - train.txt
# - validation.txt
# - test.txt
#
# =========================================================

# You should prepare the phenotype and covariates files according to your study design and analysis plan. The files should contain appropriate column names for sample IDs, phenotypes, and covariates. And the phenotype files shouls be divided into training, and testing sets according to your study design. The sample IDs in the phenotype and covariates files should match those in the genotype data for downstream analyses.

# =========================================================
# 7. GWAS summary statistics preparation
# =========================================================
# We supposed that you have downloaded the GWAS summary statistics you have,and prepare them in the format required by LDpred2
# Here we need GWAS file contains column: rsid chr pos a1 a0 beta beta_se n_eff
# And here we use snpID as chr:pos:sorted(a1:a0), which matched the SNP ID in the LD reference panel and genotype data.