# CDS-2024-Fall-Capstone


## 09/27/2024
Familiar with the backgound and go through papaers, mostly Methods.

Paper
1. https://www.nature.com/articles/s41588-024-01682-1 : Anni
2. https://www.nature.com/articles/s41588-024-01689-8 : Lyric
3. https://www.biorxiv.org/content/10.1101/2023.04.26.538501v2.full : Catherine

Past code and model:
1. SCENT (poisson regression) Github: https://github.com/immunogenomics/SCENT
2. SCAR (regression) Github: https://github.com/snehamitra/SCARlink

##  10/03/2024
Meeting to aggregate questions and information. 8 pm 

## Environment Set-Up
1. Open NYU HPC: https://ood.hpc.nyu.edu/

2. Open the HPC Console

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
rm -rf results && rm -rf resources
rm -rf *.txt && rm -rf *.out && rm -rf conda_cache && rm -rf .snakemake/conda/ && conda clean --all

# Create a specific env
conda env create -f <file_name>.yaml
conda activate <file_name>
conda deactivate

# Find top 20 files that take the most space
du -ah | sort -rh | head -n 20
```

<!-- ## Modifications

- `seurat.yaml`: Commented out `macs2`, instead loads HPC's `macs2`
- `run_snakemake.sh`: HPC has its own `snakemake` package, so no need to create one
- `run_snakemake.sh`: HPC uses a different job scheduler, so the original `bsub` command was changed to `sbatch`
- `SCENTfunctions.R`: added a `library(Matrix)` call to import package

-->