# Tutorial:
###
### 
# https://privefl.github.io/bigsnpr/articles/LDpred2.html
# https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html


# Timing: 2–4 h depending on sample size and computational resources

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

# 9 hours are requied to finished for 19 files
# Runtime depends on the trait and its GWAS characteristics, but is generally 30–60 minutes.

while read -r trait ancestry
do

    echo "Running ${trait} (${ancestry})"

    Rscript ${SCRIPT_DIR}/LDpred2.R \
        -o ${LDpred2_DIR}/${trait} \
        -g ${GWAS_DIR}/${trait}.hm3.gwas.txt.gz \
        -p ${LD_REF_DIR}/${ancestry}/${ancestry}_chr#.hm3.

done < ${GWAS_DIR}/traits.sh


# For some parameter settings, no score may be generated. In such cases, NA values are replaced with 0.
cd $LDpred2_DIR

for file in `ls *dpred2_grid.tsv.gz`
do
  zcat $file |
  sed 's/NA/0/g' |
  bgzip > tmp.gz
  mv tmp.gz $file
done


# Trouble shooting:
# h2_est < 0 ,we can estimate the h2 using LDSC instead and put the 