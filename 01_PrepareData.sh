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
J=8
PROJECT_DIR=$(pwd)

# =========================================================
# 1. Set project directories
# =========================================================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -dir)    PROJECT_DIR="$2"; shift ;;
        -j)      J="$2"; shift ;;
        -skipdl) SKIP_DL=true ;;
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

# Create directories
mkdir -p ${VCF_DIR}
mkdir -p ${GENETIC_MAP_DIR}
mkdir -p ${LD_REF_DIR}
mkdir -p ${LDpred2_DIR}
mkdir -p ${PRS_DIR}
mkdir -p ${GENO_DIR}
mkdir -p ${Layer1_DIR}
mkdir -p ${Layer2_DIR}
mkdir -p ${Weight_DIR}

# Create ancestry-specific directories
for POP in AFR EUR EAS SAS AMR
do
    mkdir -p ${LD_REF_DIR}/${POP}
done

# =========================================================
# 2. Download 1000 Genomes Phase 3 VCF files
# =========================================================
cd ${VCF_DIR}

if [[ "$SKIP_DL" != true ]]; then
    printf '%s\n' {1..22} | xargs -P ${J} -I {} wget -c -t 0 \
    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz    
fi

# Download population information file
if [[ "$SKIP_DL" != true ]]; then
    wget -c -t 0 \
    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
fi


# =========================================================
# 3. Download genetic map files
# =========================================================
#
# IMPORTANT:
#
# The directory name MUST be:
#
#   omni_genetic_map
#
# These files are required by LDpred2 for
# LD reference construction.
#
# The genetic map is based on GRCh37 coordinates.
#
# =========================================================
cd ${GENETIC_MAP_DIR}

if [[ "$SKIP_DL" != true ]]; then
    printf '%s\n' {1..22} | xargs -P ${J} -I {} wget -c -t 0 \
    https://raw.githubusercontent.com/joepickrell/1000-genomes-genetic-maps/master/interpolated_OMNI/chr{}.OMNI.interpolated_genetic_map.gz
fi
