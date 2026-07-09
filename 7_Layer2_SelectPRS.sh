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


# Get the multi-ancestry and multi-traits PRS
cd ${Layer2_DIR}

mkdir -p GPS MixPRS log

# Please change the paths in config.sh according to your own data directory structure.
echo "allpheno=\"${PHENO_DIR}/${TARGET_PHENO}.txt\"
trainpheno=\"${PHENO_DIR}/train.${TARGET_PHENO}.txt\"
testpheno=\"${PHENO_DIR}/test.${TARGET_PHENO}.txt\"
cov=\"age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10\"
cov_formula=\"age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10\"
method=\"gaussian\"
pheno=\"${TARGET_PHENO}\"" > "config.sh"

# 7-1.Normalize multi ancestry PRS
source config.sh

for trait in T2D BMI HbA1c WHR
do

    echo -e "========== ${trait} processing ============\n"

    zcat ${Layer1_DIR}/${trait}/MixPRS/${trait}-MulAnc.sscore.gz |
    KeyMapReplacer.py -k1 -a NA -p <(cat ${trainpheno}) -x |
    Rscript ${SCRIPT_DIR}/Scale.R -c ${trait}_MulAnc_PRS -t ${trait}_NormPRS |
    wcut -t "IID,${trait}_NormPRS" > GPS/${trait}-MulAnc.Norm.sscore

done


# 7-2. Perform stepAIC to select informative multi-trait PRSs.
csvtk -t join -f IID GPS/*-MulAnc.Norm.sscore |
  KeyMapReplacer.py -p<(cat $trainpheno | wcut -t "IID,${pheno}") -k1 -a NA -x |
  wcut -f2- |
  Rscript ${SCRIPT_DIR}/stepAIC.R -f "${pheno}~." -m ${method} |
  tee log/Layer2.stepAIC.log

# 7-3. Save the selected Layer 1 multi-ancestry PRSs for Layer 2 multi-trait PRS construction.

grep "NormPRS " log/Layer2.stepAIC.log |
  grep -v "+" |
  wcut -f1 |
  sed "s/_NormPRS/_MulAnc_PRS/" > stepAIC.para.sh

cat stepAIC.para.sh |
 wcut -f1 | csvtk -D"+" transpose > mixing_fomular.sh

for trait in $(sed 's/_MulAnc_PRS//g' stepAIC.para.sh)
do

    zcat ${Layer1_DIR}/${trait}/MixPRS/${trait}-MulAnc.sscore.gz > GPS/${trait}-MulAnc.raw.sscore

done

# 7-4. Fit the final Layer 2 multi-trait PRS model.

# The selected Layer 1 multi-ancestry PRSs are merged by IID and used as predictors in the final regression model.
# The mixing formula is generated from the stepAIC-selected features.
selected_scores=$(sed 's/_MulAnc_PRS//g' stepAIC.para.sh | sed 's|^|GPS/|;s|$|-MulAnc.raw.sscore|' | paste -sd' ')
n_selected=$(echo "$selected_scores" | wc -w)
formula=$(cat mixing_fomular.sh)

csvtk -t join -f IID ${selected_scores} |
  KeyMapReplacer.py -p <(cat ${trainpheno}) -k1 -a NA -x |
  Rscript ${SCRIPT_DIR}/GlmRegressiono.R -f "${pheno}~${formula}+${cov_formula}" -m ${method} -n "${pheno}~${cov_formula}" -r y |
  tee log/Layer2.Mix.log


# 7-5. Extract Layer 2 mixing coefficients and construct the final multi-ancestry, multi-trait PRS.

# Extract regression coefficients for the PRSs selected by stepAIC.
# These coefficients will be used as Layer 2 mixing weights.
cat log/Layer2.Mix.log |
  SubsetByKey.py -f<(cat stepAIC.para.sh | wcut -f1) -c1 -k |
  wcut -f1,2 > mixing.beta.sh


while read prs beta
do

    awk -v b="${beta}" '
        BEGIN{OFS="\t"}
        NR==1 {print; next}
        {print $1, $2*b}
    ' GPS/${prs%_MulAnc_PRS}-MulAnc.raw.sscore > MixPRS/${prs%_MulAnc_PRS}-mix.beta.txt

done < mixing.beta.sh

# Sum all the weighted scores into final multi-ancestry, multi-traits PRS
n=$(ls MixPRS/*-mix.beta.txt 2>/dev/null | wc -l)

if [ "$n" -eq 0 ]; then
  echo "No weighted PRS files found in MixPRS/"
  exit 1
fi
if [ "$n" -gt 1 ]; then
  csvtk -t join -f IID MixPRS/*-mix.beta.txt
else
  cat MixPRS/*-mix.beta.txt

fi |
awk '
BEGIN{OFS="\t"}
NR==1 {next}
{
  prs=0
  for(i=2;i<=NF;i++) prs += $i
  print $1, prs
}
' |
AddTitle.py "IID PRS" |
csvtk space2tab |
bgzip > MixPRS/${pheno}-MulAnc.MultiTraits.sscore.gz