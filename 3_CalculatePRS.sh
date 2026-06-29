# Calculate PRS
PROJECT_DIR=$(pwd)

# =========================================================
# 1. Set project directories
# =========================================================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -dir)    PROJECT_DIR="$2"; shift ;;
        -skip)    SKIPCMP=true; shift ;;
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

cd $PRS_DIR

if [[ "$SKIPCMP" != true ]]; then
    for trait in $(cut -f1 ${GWAS_DIR}/traits.sh)
    do
        for chr in {1..22}
        do

            plink2 \
                --pfile ${GENO_DIR}/EUR-${chr} \
                --score <(zcat ${LDpred2_DIR}/${trait}ldpred2_grid.tsv.gz | body grep "^${chr}:") 1 2 header-read zs list-variants-zs cols='sid,nallele,dosagesum,scoresums' \
                --score-col-nums 4-105 \
                --memory 4096 \
                --threads 8 \
                --out ${PRS_DIR}/${trait}-chr${chr}

        done
    done
fi


# Sum scores for all chromosomes
for trait in $(cut -f1 ${GWAS_DIR}/traits.sh)
do

    zstdcat ${PRS_DIR}/${trait}-chr*.sscore.zst |
        RemoveDuplicateTitle.py |
        wcut -f1,5- |
        body datamash -s -g 1 sum 2-103 |
        sed 's/#IID/IID/' |
        bgzip > ${PRS_DIR}/${trait}-sscore.gz

done

# rm ${PRS_DIR}/*-chr*

# Save model parameters for downstream analyses.
Rscript -e '
p_seq <- signif(exp(seq(log(1e-4), log(1), length.out=17)), 2)

for(sparse in c("F","T")){
  for(h2 in c(0.7,1,1.4)){
    for(p in p_seq){
      cat(sprintf("p:%s_h2:%s_%s_SUM\n", p, h2, sparse))
    }
  }
}
' > parameters.sh