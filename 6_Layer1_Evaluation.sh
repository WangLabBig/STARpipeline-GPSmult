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

# 6-1. Evaluate the performance of the multi-ancestry PRS in the training set.
# perform this in Layer1/xxxpheno dir
# Steps: regress out PCs -> normalize PRS -> fit regression model.
for tt in BMI T2D WHR HbA1c
do
    cd ${Layer1_DIR}/${tt}

    source config.sh
    prs=$(cat PRS.sh)

    zcat MixPRS/${prs}-MulAnc.sscore.gz |
        KeyMapReplacer.py -k1 -a NA -p<(cat ${trainpheno}) -x |
        wcut -t "${pheno},${prs}_MulAnc_PRS,${cov}" |
        Rscript ${SCRIPT_DIR}/Residuals.R -f "${prs}_MulAnc_PRS~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10" -t adjPRS |
        Rscript ${SCRIPT_DIR}/Scale.R -c adjPRS -t adjNormPRS |
        Rscript ${SCRIPT_DIR}/GlmRegression.R -f "${pheno}~adjNormPRS+${cov_formula}" -m ${method} -n "${pheno}~${cov_formula}" -r y |
        tee log/${prs}-multi.train.log

    # 6-2. Evaluate the performance of the multi-ancestry PRS in the testing set.
    zcat MixPRS/${prs}-MulAnc.sscore.gz |
        KeyMapReplacer.py -k1 -a NA -p<(cat ${testpheno}) -x |
        wcut -t "${pheno},${prs}_MulAnc_PRS,${cov}" |
        Rscript ${SCRIPT_DIR}/Residuals.R -f "${prs}_MulAnc_PRS~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10" -t adjPRS |
        Rscript ${SCRIPT_DIR}/Scale.R -c adjPRS -t adjNormPRS |
        Rscript ${SCRIPT_DIR}/GlmRegression.R -f "${pheno}~adjNormPRS+${cov_formula}" -m ${method} -n "${pheno}~${cov_formula}" -r y |
        tee log/${prs}-multi.test.log

done
# ===================== End of Layer 1 =====================

# At this stage, the multi-ancestry PRS has been trained and evaluated for a single target trait.

# Repeat Steps 4–6 for all remaining traits in the Layer1 directory.

# After completing Layer 1 for all traits, proceed to Layer 2 to construct and evaluate the multi-trait PRS model.