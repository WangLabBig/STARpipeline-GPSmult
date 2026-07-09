# Use the optimal adjusted and normalized PRSs (adjNormPRS) from Layer 1.
# Repeat the stepAIC procedure for all traits in the Layer 1 directory.

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


# 5-1. Perform feature selection using stepAIC to identify informative PRSs.
for tt in BMI T2D WHR HbA1c
do
    cd ${Layer1_DIR}/${tt}

    source config.sh

    prs=$(cat PRS.sh)
    trait_list=$(paste -sd, traits.sh)
    n=$(ls GPS/*-adjNorm.best.score.txt 2>/dev/null | wc -l)

    if [ "$n" -gt 1 ]; then

        csvtk -t join -f IID GPS/*-adjNorm.best.score.txt |
        KeyMapReplacer.py -p<(cat ${trainpheno}) -k1 -a NA -x |
        wcut -t "${pheno},${trait_list}" |
        Rscript ${SCRIPT_DIR}/stepAIC.R -f "${pheno}~." -m ${method} |
        tee log/${prs}.stepAIC.log

    elif [ "$n" -eq 1 ]; then

        cat GPS/*-adjNorm.best.score.txt |
        KeyMapReplacer.py -p<(cat ${trainpheno}) -k1 -a NA -x |
        wcut -t "${pheno},${trait_list}" |
        Rscript ${SCRIPT_DIR}/stepAIC.R -f "${pheno}~." -m ${method} |
        tee log/${prs}.stepAIC.log

    else

        echo "No *-adjNorm.best.score.txt files found in GPS/"

    fi


    # 5-2. Save parameters that survival after stepAIC
    grep "^${prs}" log/${prs}.stepAIC.log |
    wcut -f1 |
    KeyMapReplacer.py -p<(cat best_beta_para.sh) -a NA -k1 -x |
    wcut -f1,2 > stepAIC.para.sh


    # 5-3. Generate regression fomular for determining mixing weights
    cat stepAIC.para.sh |
    wcut -f1 | csvtk -D"+" transpose > mixing_fomular.sh


    # 5-4. Fit the final Layer 1 mixing model using the features selected by stepAIC.
    # If multiple PRSs are selected, they are merged by IID before model fitting.
    # If only one PRS is selected, the score file is used directly.
    formula=$(cat mixing_fomular.sh)

    selected_scores=$(cat stepAIC.para.sh | wcut -f1 | sed 's|^|GPS/|;s|$|-raw.best.score.txt|' | paste -sd' ')
    n_selected=$(echo "$selected_scores" | wc -w)

    if [ "$n_selected" -gt 1 ]; then

        csvtk -t join -f IID ${selected_scores} |
        KeyMapReplacer.py -p<(cat ${trainpheno}) -k1 -a NA -x |
        Rscript ${SCRIPT_DIR}/GlmRegressiono.R -f "${pheno}~${formula}+${cov_formula}" -m ${method} -n "${pheno}~${cov_formula}" -r y |
        tee log/${prs}.layer1Mix.log

    elif [ "$n_selected" -eq 1 ]; then

        cat ${selected_scores} |
        KeyMapReplacer.py -p<(cat ${trainpheno}) -k1 -a NA -x |
        Rscript ${SCRIPT_DIR}/GlmRegressiono.R -f "${pheno}~${formula}+${cov_formula}" -m ${method} -n "${pheno}~${cov_formula}" -r y |
        tee log/${prs}.layer1Mix.log

    else

        echo "No selected scores found in stepAIC.para.sh"

    fi


    # 5.5 Extract mixing beta from multiple ancestry linear regression.
    cat log/${prs}.layer1Mix.log |
    SubsetByKey.py -f<(cat stepAIC.para.sh | wcut -f1) -c1 -k |
    wcut -f1,2 > mixing.beta.sh

    # 5-6. Construct the multi-ancestry PRS.

    # Multiply each ancestry-specific PRS by its corresponding mixing coefficient, estimated from the Layer 1 mixing model.
    # Sum all weighted PRSs to generate the final multi-ancestry PRS.
    # The resulting weighted PRSs are saved in the MixPRS directory.

    while read trait beta
    do

        awk -v b="$beta" '
            BEGIN{OFS="\t"}
            NR==1 {print $1,$2; next}
            {print $1,$2*b}
        ' GPS/${trait}-raw.best.score.txt > MixPRS/${trait}-mix.beta.txt

    done < mixing.beta.sh


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
    AddTitle.py "IID ${prs}_MulAnc_PRS" |
    csvtk space2tab |
    bgzip > MixPRS/${prs}-MulAnc.sscore.gz
done
