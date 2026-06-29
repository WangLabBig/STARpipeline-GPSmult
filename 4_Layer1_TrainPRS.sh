# Layer 1 consists of three steps:
# (1) Select the optimal parameter set for each ancestry-specific PRS.
# (2) Perform feature selection to retain informative PRSs.
# (3) Linearly combine the selected PRSs.
# This process finally generates multi-ancestry PRS foa a specific trait.
PROJECT_DIR=$(pwd)
TARGET_PHENO="BMI"
# =========================================================
# 1. Set project directories
# =========================================================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -dir)    PROJECT_DIR="$2"; shift ;;
        -p)      TARGET_PHENO="$2"; shift ;;
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


cd ${Layer1_DIR}

# Create directories for each trait.
# The traits listed below are examples; please customize them according to your analysis.
mkdir -p BMI HbA1c T2D WHR


# Each trait directory should include a config.sh file with the following contents:
# Update the paths in config.sh to reflect your own directory structure and file locations.
for i in BMI HbA1c T2D WHR
do
  echo "PRS_dir=\"${PRS_DIR}\"
allpheno=\"${PHENO_DIR}/${TARGET_PHENO}.txt\"
trainpheno=\"${PHENO_DIR}/train.${TARGET_PHENO}.txt\"
testpheno=\"${PHENO_DIR}/test.${TARGET_PHENO}.txt\"
cov=\"age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10\"
cov_formula=\"age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10\"
method=\"gaussian\"
pheno=\"${TARGET_PHENO}\"" > "${i}/config.sh"

  mkdir -p ${i}/train_result;
  mkdir -p ${i}/MixPRS;
  mkdir -p ${i}/log;
  mkdir -p ${i}/GPS;
    
  echo ${i} > ${i}/PRS.sh;
  cp ${PRS_DIR}/parameters.sh ${i}/parameters.sh;
  grep ${i} ${GWAS_DIR}/traits.sh | cut -f1 > ${i}/traits.sh;
done



# 4-1.Run training for each trait and each parameter:
# repeat the training step for all your traits in the layer1 directory, and save the training results in the train_result directory under each trait's directory. You can also change the method (gaussian, binomial) in the config.sh file for each trait according to your phenotype type.
for tt in BMI T2D WHR HbA1c
do
    cd ${Layer1_DIR}/${tt}
    source ${Layer1_DIR}/${tt}/config.sh
    while read trait
    do

        while read parameter
        do

            zcat ${PRS_DIR}/${trait}-sscore.gz |

            KeyMapReplacer.py -k1 -a NA -p <(cat ${trainpheno}) -x |

            sed "s/${parameter}/PRS/g" |

            # regress out top 10 PCs
            Rscript ${SCRIPT_DIR}/Residuals.R -f 'PRS~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10' -t adjPRS |

            # normalize adjPRS
            Rscript ${SCRIPT_DIR}/Scale.R -c adjPRS -t adjNormPRS |

            wcut -t "${pheno},adjNormPRS,${cov}" |
            Rscript ${SCRIPT_DIR}/GlmRegression.R -f "${pheno}~adjNormPRS+${cov_formula}" -m ${method} -n "${pheno}~${cov_formula}" -r y > train_result/${trait}-${parameter}.train.txt

        done < parameters.sh

    done < traits.sh



# 4-2.Select the best performing score for cohort-specific, ancestry-stratified PRS.
# Also repeat this step for all your traits in the layer1 directory
    for trait in $(cat traits.sh)
    do

        grep 'Incremental_Model_Rsq' train_result/${trait}-*.txt |
        sed 's/:Incremental_Model_Rsq://' |
        sed 's|train_result/||' |
        bgzip > train_result/${trait}-train.R2.summary.txt.gz;

        zcat train_result/${trait}-train.R2.summary.txt.gz |
        wcut -f1,2 |
        datamash -f max 2 |
        wcut -f1,2 > train_result/${trait}-train.best.R2.txt

    done


# 4-3. Identify the optimal PRS parameter set for each ancestry based on the maximum R².
    cat train_result/*best*R2*txt |
    sed 's/-/ /' |
    sed 's/.train.txt//' |
    csvtk space2tab > best_R2_para.sh


# 4-4. Extract the corresponding effect size (beta coefficient) for the selected PRS.
    while read trait para _
    do

        grep 'adjNormPRS ' train_result/${trait}-${para}.train.txt |
        sed "s/adjNormPRS/${trait} ${para}/" |
        wcut -f1-3 |
        csvtk space2tab >> best_beta_para.sh

    done < best_R2_para.sh


# 4-5. Generate the optimal PRS for each trait and ancestry, and save the scores in the GPS directory.

# Two versions of the PRS are generated:
# (1) the raw PRS without adjustment, and
# (2) the adjusted and normalized PRS (adjNormPRS) used for model training.

# The raw PRS is retained to enable reconstruction of SNP-level weights using the raw score mixing weights.

    while read trait para _
    do

        # Generate adjusted and normalized PRS for the training set
        zcat ${PRS_DIR}/${trait}-sscore.gz |
        KeyMapReplacer.py -p <(cat ${trainpheno}) -k1 -a NA -x |
        sed "s/${para}/PRS/g" |
        Rscript ${SCRIPT_DIR}/Residuals.R -f 'PRS~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10' -t adjPRS |
        Rscript ${SCRIPT_DIR}/Scale.R -c adjPRS -t adjNormPRS |
        sed "s/adjNormPRS/${trait}/g" |
        wcut -t "IID,${trait}" > GPS/${trait}-adjNorm.best.score.txt

        # Generate raw PRS for all individuals
        zcat ${PRS_DIR}/${trait}-sscore.gz |
        KeyMapReplacer.py -p <(cat ${allpheno}) -k1 -a NA -x |
        sed "s/${para}/${trait}/g" |
        wcut -t "IID,${trait}" > GPS/${trait}-raw.best.score.txt

    done < best_beta_para.sh

done

