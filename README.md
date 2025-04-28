# Machine Learning for Multiomics: Identifying Interacting DNA Elements

Authors: Anni Zheng, Catherine Hua, Lyric Li  
Affiliation: NYU Center for Data Science  
Course: DS-1006 Capstone, Fall 2024  

## Project Overview

Enhancers play a pivotal role in regulating gene expression and often interact with each other to modulate transcriptional activity. However, the mechanisms and prevalence of enhancer–enhancer interactions remain poorly understood.

In this project, we developed a scalable, fully functional pipeline under the NYU High Performance Computing (HPC) environment to systematically identify and analyze potential enhancer interactions:
	•	We integrated a Generalized Linear Model (GLM) and Poisson regression framework.
	•	We analyzed enhancer accessibility and gene expression at the single-cell level.
	•	We focused our analysis on the CD14-Mono cell type from a PBMC dataset.

Our results suggest that enhancer activities may exhibit non-linear, saturated behaviors when acting in combination, highlighting the complexity of enhancer interactions.

## Dataset

The input data includes:
	•	atac_matrix: Chromatin accessibility (0–1 scale) for enhancer regions across cells.
	•	rna_matrix: Single-cell RNA-seq gene expression counts (0–5585).
	•	meta_data: Metadata including cell type (CD14-Mono), total UMIs, and mitochondrial content percentage.

Data Preprocessing:
	•	Filtered genes and enhancers retained only if >5% of cells had nonzero counts.
	•	Focused analysis on CD14-Mono cell type.


## Methodology

SCENT Model
	•	Poisson regression to associate enhancer accessibility with gene expression.
	•	Applied bootstrap-based significance testing and FDR correction (threshold < 0.1).

Epistasis Model
	•	Built three models:
	•	cells10: Enhancer 1 active, Enhancer 2 inactive.
	•	cells01: Enhancer 1 inactive, Enhancer 2 active.
	•	cells11: Both enhancers active.
	•	Compared models against the baseline cells00 (both enhancers inactive).
	•	Evaluated interaction strength via beta coefficients and statistical significance.

## Pipeline Workflow
	1.	Gene and Peak Annotation
	•	Download GENCODE annotations.
	•	Annotate genes and enhancer peaks.
	•	Match genes and peaks within 500kb window using bedtools.
	2.	Gene and Peak Filtering
	•	Select highly variable genes (via Seurat’s MVP method).
	•	Split enhancer–gene pairs for parallelized processing.
	3.	Run SCENT Analysis
	•	Identify significant enhancer–gene pairs.
	4.	Run Epistasis Model
	•	Model enhancer–enhancer interactions.
	•	Bootstrap results and apply multiple testing corrections.

HPC Setup:
	•	Executed on NYU Greene cluster using Slurm scheduler.
	•	Snakemake workflow engine for orchestration.
	•	Used custom conda environments to handle dependency issues.

## Results and Key Findings
| Analysis                      | Key Insight |
|--------------------------------|-------------|
| SCENT output QQ Plot           | Evidence of significant enhancer–gene associations. |
| Volcano plots (Epistasis)      | Positive beta estimates stronger when both enhancers are active. |
| Beta Sum vs Beta Both          | Summed single enhancer effects often exceed dual enhancer effects — suggests **sub-additive interactions** and possible **regulatory saturation**. |

## Limitations
	•	48-hour runtime limit on NYU HPC constrained large-scale analyses.
	•	No GPU acceleration for GLM and SCENT models.
	•	Instability issues occasionally encountered with Snakemake’s DAG planning.
	•	Current analysis focused on one cell type (CD14-Mono).

## Future Work
	•	Expand analysis to all 30 PBMC cell types.
	•	Explore non-linear models (e.g., deep learning, tree-based methods) for complex interaction patterns.
	•	Investigate synthetic augmentation of scRNA-seq and scATAC-seq data to better handle data sparsity.
 
## References
1. https://www.nature.com/articles/s41588-024-01682-1 
2. https://www.nature.com/articles/s41588-024-01689-8
3. https://www.biorxiv.org/content/10.1101/2023.04.26.538501v2.full

Past codes and models:
1. SCENT (poisson regression) Github: https://github.com/immunogenomics/SCENT
2. SCAR (regression) Github: https://github.com/snehamitra/SCARlinka

## Environment Set-Up
1. NYU HPC: https://ood.hpc.nyu.edu/
2. Open HPC Console
3. Clone or Update the Repository

Clone
```
git clone https://github.com/ANNIZHENG/CDS-2024-Fall-Capstone.git
```

Update
```
git pull origin main
```

4. Manually create data folder and upload data

The path should match what's specified in CDS-2024-Fall-Capstone/config/config_PBMC.yaml

```
cd CDS-2024-Fall-Capstone
mkdir data
```

5. Send task to HPC Cluster

Direct Cluster Training
```
sbatch --time=48:00:00 run_snakemake.sh
```

Interactive Session
```
srun -t 1:00:00 --mem 100G --pty /bin/bash
bash run_snakemake.sh
```

5. Useful command to clear some memory space

```
# Remove all created files from the pipeline
rm -rf resources && rm -rf results 
rm -rf *.txt && rm -rf *.out && rm -rf *.err
rm -rf conda_cache && rm -rf .snakemake/conda/ && conda clean --all

# Create a specific env
conda env create -f <file_name>.yaml
conda activate <file_name>
conda deactivate

# Find top 20 files that take the most space
du -ah | sort -rh | head -n 20
```

<!-- ## Modifications

- `seurat.yaml`: Commented out `macs2`, instead loads HPC's `macs2`
- `SCENTfunctions.R`: added a `library(Matrix)` call to import package

-->
