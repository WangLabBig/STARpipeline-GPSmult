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

# 9-1. Extract the layer1 and layer 2 mixing weight for each score file:
# Define your layer 1 directory:
cd $Weight_DIR

mkdir -p mixing_weights


# Target trait:
target_trait=BMI

for i in BMI T2D HbA1c WHR
do
    cat  ${Layer1_DIR}/${i}/best_R2_para.sh >> R2_performance.txt;
    cat  ${Layer1_DIR}/${i}/best_beta_para.sh >> beta_performance.txt;
    cat ${Layer1_DIR}/${i}/mixing.beta.sh >> layer1.mix.weight.txt;
done

# Layer 1 mixing weights
cat ${Layer2_DIR}/mixing.beta.sh |
  sed 's/_MulAnc_PRS//g' > layer2.mixing.beta.txt;

cat beta_performance.txt |
  KeyMapReplacer.py -p R2_performance.txt -k1,2 -a NA -x |
  sed 's/_SUM//g' |
  AddTitle.py 'PRS BestParameter Beta deltaR2' |
  KeyMapReplacer.py -p<(cat layer1.mix.weight.txt | AddTitle.py 'PRS Layer1.mix.beta') -k1 -a NA > layer1.mix.res.txt


# Layer 2 mixing weights
cat layer1.mix.res.txt |
  awk '{split($1,a,".");print a[1],$0}' |
  KeyMapReplacer.py -p <(cat layer2.mixing.beta.txt | AddTitle.py 'PRS Layer2.mix.beta') -k1 -a NA |
  sed 's/PRS/Trait/' |
  sed 's/NA/0/g' |
  awk 'BEGIN{OFS="\t"} NR==1{
    print $0,"FinalWeights"
    next}
    {print $0,$6*$7}' >  mixing_weights/${target_trait}.PRS.mixing.weight.txt

rm *txt

cat mixing_weights/${target_trait}.PRS.mixing.weight.txt |
    wcut -t 'PRS,BestParameter,FinalWeights' | 
    awk '$3!=0' |
    tail -n+2  |
    cut -f1,2 > PRS.para.txt


cat mixing_weights/${target_trait}.PRS.mixing.weight.txt |
    wcut -t 'PRS,BestParameter,FinalWeights' | 
    awk '$3!=0' |
    tail -n+2 |
    cut -f1,3 > PRS.weights.txt


# 9-2. Extract the snp weight file.
# Only need to use the file with final mixing weights not equals 0.  
LDpred2_DIR=${DATA_DIR}/LDpred2
hm3=${HM3_DIR}/hapmap3.snp

while read prs para
do

  zcat ${LDpred2_DIR}/${prs}ldpred2_grid.tsv.gz |
  wcut -t "rsID,a1,${para}" |
  tail -n +2 |
  AddTitle.py "snpID a1 beta" |
  csvtk space2tab > ${prs}.txt

done < PRS.para.txt


cat $hm3 |
  awk 'BEGIN {FS="\t";OFS="\t"} {split($1, a, ":");print $0,a[4]}' |  
  AddTitle.py 'rsID EA' |
  csvtk space2tab > snpid.txt


Rscript ${SCRIPT_DIR}/SNP.weights.R \
    snpid.txt \
    PRS.weights.txt \
    ${target_trait}.PRS.weight.txt 


# ==============================================================================


############# Test the weight calculation is correct or not ##############

mkdir -p GPS

for chr in {1..22}
do

  plink2 \
    --pfile ${GENO_DIR}/EUR-${chr} \
    --score <(cat ${target_trait}.PRS.weight.txt | body grep "^${chr}:") 1 2 header-read zs list-variants-zs cols='sid,nallele,dosagesum,scoresums' \
    --score-col-nums 3 \
    --memory 4096 \
    --threads 8 \
    --out GPS/${target_trait}-chr${chr}

done


# Sum scores for all chromosomes
# parallel -j1 -k "
# zstdcat GPS/{1}-chr*.sscore.zst |
#     RemoveDuplicateTitle.py |
#     wcut -f1,5 |
#     body datamash -s -g 1 sum 2 |
#     sed 's/#IID/IID/' |
#     bgzip > GPS/{1}-sscore.gz
# " ::: `echo ${target_trait}`


zstdcat GPS/${target_trait}-chr*.sscore.zst |
        RemoveDuplicateTitle.py |
        wcut -f1,5 |
        body datamash -s -g 1 sum 2 |
        sed 's/#IID/IID/' |
        bgzip > GPS/${target_trait}-sscore.gz

## Calculate the Pearson correlation between beta-merged score and PRS mixing score.

mix_PRS="${PROJECT_DIR}/Layer2/MixPRS/BMI-MulAnc.MultiTraits.sscore.gz"

echo 'the pearson R between beta-merged score and PRSs mixing score:'
zcat GPS/${target_trait}-sscore.gz |
  KeyMapReplacer.py -p<(zcat ${mix_PRS}) -k1 -a NA -x |
  Rscript ${SCRIPT_DIR}/PearsonR2.R -f final_beta_SUM -s PRS 