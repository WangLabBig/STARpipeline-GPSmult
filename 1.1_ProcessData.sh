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
DATA_DIR=${PROJECT_DIR}/data
GWAS_DIR=${DATA_DIR}/GWAS
HM3_DIR=${DATA_DIR}/HapMap3
hm3=${HM3_DIR}/hapmap3.snp
GENO_DIR=${DATA_DIR}/Genotype
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
# 1. Generate ancestry-specific sample lists
# =========================================================
cd ${VCF_DIR}
PANEL=${VCF_DIR}/integrated_call_samples_v3.20130502.ALL.panel

for POP in AFR EUR EAS SAS AMR
do
    awk -v p=${POP} '$3==p {print $1, $1}' ${PANEL} > ${LD_REF_DIR}/${POP}/${POP}.keep
done


# =========================================================
# 2. Convert VCF files to PLINK format
# =========================================================
#
# Output:
#   ${LD_REF_DIR}/${POP}/${POP}_chr${CHR}
#
# =========================================================
for POP in AFR EUR EAS SAS AMR
do
    for CHR in {1..22}
    do
        plink2 \
            --vcf ${VCF_DIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
            --keep ${LD_REF_DIR}/${POP}/${POP}.keep \
            --double-id \
            --max-alleles 2 \
            --maf 0.05 \
            --make-bed \
            --out ${LD_REF_DIR}/${POP}/${POP}_chr${CHR}
    done
done


# =========================================================
# 3. Recode SNP IDs
# =========================================================
#
# SNP IDs are reformatted as:
#
#   CHR:POS:A1:A2
#
# Alleles are alphabetically sorted to ensure
# consistent SNP naming across datasets.
#
# =========================================================
for POP in AFR EUR EAS SAS AMR
do
    for CHR in {1..22}
    do
        BIM=${LD_REF_DIR}/${POP}/${POP}_chr${CHR}.bim
        awk '{
            if ($5 < $6) {
                print $1, $1":"$4":"$5":"$6, $3, $4, $5, $6
            } else {
                print $1, $1":"$4":"$6":"$5, $3, $4, $5, $6
            }
        }' ${BIM} > ${BIM}.tmp
        mv ${BIM}.tmp ${BIM}
    done
done


# =========================================================
# 4. Extract HapMap3 variants
# =========================================================
for POP in AFR EUR EAS SAS AMR
do
    for CHR in {1..22}
    do
        plink2 \
            --bfile ${LD_REF_DIR}/${POP}/${POP}_chr${CHR} \
            --extract ${hm3} \
            --make-bed \
            --out ${LD_REF_DIR}/${POP}/${POP}_chr${CHR}.hm3
    done
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