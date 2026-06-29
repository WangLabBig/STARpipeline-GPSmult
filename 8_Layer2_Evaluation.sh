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


# 8-1. Evaluate the performance of multi-ancestry, multi-traits PRS in the training set.
# Steps: regress out PCs -> normalize PRS -> fit regression model.
cd ${Layer2_DIR}
source config.sh

zcat MixPRS/${pheno}-MulAnc.MultiTraits.sscore.gz |
  KeyMapReplacer.py -k1 -a NA -p<(cat $trainpheno) -x |
  wcut -t "${pheno},PRS,$cov" |
  Rscript ${SCRIPT_DIR}/Residuals.R -f 'PRS~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10' -t  adjPRS |
  Rscript ${SCRIPT_DIR}/Scale.R -c adjPRS -t adjNormPRS |
  Rscript ${SCRIPT_DIR}/GlmRegression.R -f "${pheno}~adjNormPRS+${cov_formula}" -m ${method} -n "${pheno}~${cov_formula}" -r y |
  tee log/Layer2.train.log


# 8-2. Evaluate multi-ancestry mixed PRS performance on the testing set.
# Steps: regress out PCs -> normalize PRS -> fit regression model.
zcat MixPRS/${pheno}-MulAnc.MultiTraits.sscore.gz |
  KeyMapReplacer.py -k1 -a NA -p<(cat $testpheno) -x |
  wcut -t "${pheno},PRS,$cov" |
  Rscript ${SCRIPT_DIR}/Residuals.R -f 'PRS~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10' -t  adjPRS |
  Rscript ${SCRIPT_DIR}/Scale.R -c adjPRS -t adjNormPRS |
  Rscript ${SCRIPT_DIR}/GlmRegression.R -f "${pheno}~adjNormPRS+${cov_formula}" -m ${method} -n "${pheno}~${cov_formula}" -r y |
  tee log/Layer2.test.log