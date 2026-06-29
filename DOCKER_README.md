# Docker Setup for PRS Pipeline

## Prerequisites

Your repository should have this layout before building:

```
project/
├── Dockerfile
├── scripts/              ← custom utility scripts from 00_Tools.sh
│   ├── wcut
│   ├── SubsetByKey.py
│   ├── KeyMapReplacer.py
│   ├── AddColumn.py
│   ├── AddTitle.py
│   ├── body
│   ├── RemoveDuplicateTitle.py
│   ├── LDpred2.R
│   ├── Ldpred2LDRefCache.R
│   ├── GlmRegression.R
│   ├── Residuals.R
│   ├── Scale.R
│   ├── stepAIC.R
│   ├── PearsonR2.R
│   └── SNP.weights.R
├── 00_Tools.sh
├── 01_PrepareData.sh
├── 1_ProcessData.sh
├── 2_RunLDpred2.sh
├── 3_CalculatePRS.sh
├── 4_Layer1_TrainPRS.sh
├── 5_Layer1_SelectPRS.sh
├── 6_Layer1_Evaluation.sh
├── 7_Layer2_SelectPRS.sh
├── 8_Layer2_Evaluation.sh
└── 9_TraceMixingWeights.sh
```

---

## Build the image

```bash
docker build -t prs-pipeline:latest .
```

Build will take ~15–30 min mainly due to R package compilation (bigsnpr).

---

## Run the pipeline

Mount your data directory into `/data` inside the container:

```bash
docker run --rm -it \
    -v /your/local/data:/data \
    prs-pipeline:latest \
    bash /opt/pipeline/01_PrepareData.sh -dir /data -j 8
```

### Run step by step

```bash
# Step 1 – Prepare data & LD reference
docker run --rm -it -v /your/local/data:/data prs-pipeline:latest \
    bash /opt/pipeline/01_PrepareData.sh -dir /data

# Step 2 – Run LDpred2
docker run --rm -it -v /your/local/data:/data prs-pipeline:latest \
    bash /opt/pipeline/2_RunLDpred2.sh -dir /data

# Step 3 – Calculate PRS
docker run --rm -it -v /your/local/data:/data prs-pipeline:latest \
    bash /opt/pipeline/3_CalculatePRS.sh -dir /data

# Step 4 – Layer 1: train PRS
docker run --rm -it -v /your/local/data:/data prs-pipeline:latest \
    bash /opt/pipeline/4_Layer1_TrainPRS.sh -dir /data

# Step 5 – Layer 1: feature selection
docker run --rm -it -v /your/local/data:/data prs-pipeline:latest \
    bash /opt/pipeline/5_Layer1_SelectPRS.sh -dir /data

# Step 6 – Layer 1: evaluate
docker run --rm -it -v /your/local/data:/data prs-pipeline:latest \
    bash /opt/pipeline/6_Layer1_Evaluation.sh -dir /data

# Step 7 – Layer 2: feature selection
docker run --rm -it -v /your/local/data:/data prs-pipeline:latest \
    bash /opt/pipeline/7_Layer2_SelectPRS.sh -dir /data

# Step 8 – Layer 2: evaluate
docker run --rm -it -v /your/local/data:/data prs-pipeline:latest \
    bash /opt/pipeline/8_Layer2_Evaluation.sh -dir /data

# Step 9 – Trace mixing weights
docker run --rm -it -v /your/local/data:/data prs-pipeline:latest \
    bash /opt/pipeline/9_TraceMixingWeights.sh -dir /data
```

---

## Notes

- **PLINK2 version**: The Dockerfile pins the Jan 2024 alpha5 build. Update the URL in step 3 if needed.
- **R version**: Installs from CRAN40 (R ≥ 4.3) on Ubuntu 22.04.
- **bigsnpr** requires the genetic map directory (`omni_genetic_map`) to be present in the working directory at runtime — mount it via `-v`.
- Hard-coded paths in your `.sh` scripts (e.g. `/hwmaster/chenqiuli/...`) should be replaced with `-dir /data` arguments or updated in each script's variable block.
- For HPC/Slurm environments, consider building a Singularity image from this Docker image:
  ```bash
  singularity pull prs-pipeline.sif docker://prs-pipeline:latest
  ```
