If you use this protocol in your research work, please cite the following papers:
- Kim, M.S., Chen, Q., Sui, Y., Yang, X., Wang, S., Weng, L.-C., Cho, S.M.J., Koyama, S., Zhu, X., Yu, K., et al. (2026). Metabolic polygenic risk scores for prediction of obesity, type 2 diabetes, and related morbidities. Cell Metab. https://doi.org/10.1016/j.cmet.2026.02.009.
- Patel, A.P., Wang, M., Ruan, Y., Koyama, S., Clarke, S.L., Yang, X., Tcheandjieu, C., Agrawal, S., Fahed, A.C., Ellinor, P.T., et al. (2023). A multi-ancestry polygenic risk score improves risk prediction for coronary artery disease. Nat. Med. 29, 1793–1803. https://doi.org/10/gsfv48.

Example data, code, and offline Docker image archive (.tar) and a Singularity image (.sif) are available for download from Zenodo (https://doi.org/10.5281/zenodo.20841307).

=====================================================================
### Prepare

**Established docker/singularity images**:
- docker:
- singularity:

**Data download**
- Example Genotype, Phenotype, Hapmap3, GWAS can be downloaded from:
- Reference data see below: step 01.download data

======================================================================
### Usage

**Run pipeline: singularity as example**
```bash
PROJ_DIR=""
SIF_DIR=""

mkdir -p ${PROJ_DIR}
mkdir -p ${PROJ_DIR}/data/
cd ${PROJ_DIR}/data/

# !!! Download example Genotype, Phenotype, Hapmap3, GWAS here
# Make sure your ${PROJ_DIR}/data/ has 4 folders:
# Genotype/  GWAS/  Phenotype/  HapMap3/

# 01.download data
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/01_PrepareData.sh -dir /data

# 1.1 process data
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/1.1_ProcessData.sh -dir /data

# 1.2 Generate LD
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/1.2_GenerateLD.sh -dir /data

# 2. Run LDpred2
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/2_RunLDpred2.sh -dir /data

# 3. Calculate PRS
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/3_CalculatePRS.sh -dir /data

# 4. Layer1 Train
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/4_Layer1_TrainPRS.sh -dir /data

# 5. Layer1 select PRS
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/5_Layer1_SelectPRS.sh -dir /data

# 6. Layer1 evaluation
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/6_Layer1_Evaluation.sh -dir /data

# 7. Layer2 Select PRS
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/7_Layer2_SelectPRS.sh -dir /data

# 8. Layer2 Evaluation
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/8_Layer2_Evaluation.sh -dir /data

# 9. Final weights
singularity exec --bind ${PROJ_DIR}:/data ${SIF_DIR}/gpsmult.sif bash /opt/pipeline/9_TraceMixingWeights.sh -dir /data
```

======================================================================

**Run pipeline: docker as example**
```bash
# load docker image: download gpsmult.tar from link above
docker load -i gpsmult.tar

# check
docker images
```

**Run analysis**
```bash
PROJ_DIR=""
DOCKER_IMAGE="gpsmult:latest"

mkdir -p ${PROJ_DIR}
mkdir -p ${PROJ_DIR}/data/
cd ${PROJ_DIR}/data/

# Make sure your ${PROJ_DIR}/data/ has 4 folders:
# Genotype/  GWAS/  Phenotype/  HapMap3/

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/01_PrepareData.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/1.1_ProcessData.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/1.2_GenerateLD.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/2_RunLDpred2.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/3_CalculatePRS.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/4_Layer1_TrainPRS.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/5_Layer1_SelectPRS.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/6_Layer1_Evaluation.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/7_Layer2_SelectPRS.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/8_Layer2_Evaluation.sh -dir /data

docker run --rm -v ${PROJ_DIR}:/data ${DOCKER_IMAGE} \
  bash /opt/pipeline/9_TraceMixingWeights.sh -dir /data
```
