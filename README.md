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

##  10/03/2004
Meeting to aggregate questions and information. 8 pm 

## Environment Set-Up
1. Open NYU HPC: https://ood.hpc.nyu.edu/

2. Open the HPC Console through **Scratch** 

Not Home* - In the website, You would find the *Files* tab in the top left, then click *Scratch*, then open the console from there.

3. Clone or Update the Repository

Clone
```
git clone https://github.com/ANNIZHENG/CDS-2024-Fall-Capstone.git
```

Update
```
git pull origin main
```

4. Send task to HPC Cluster

Direct Cluster Training
```
sbatch --time=04:00:00 run_snakemake.sh
```

Interactive Session
```
srun -t 4:00:00 --mem=100G --pty /bin/bash
bash run_snakemake.sh
```

5. Useful command to clear some memory space

```
rm -rf results && rm -rf resources
rm -rf *.txt && rm -rf *.out && rm -rf conda_cache && rm -rf .snakemake/conda/ && conda clean --all
```

<!-- ## Modifications

- `seurat.yaml`: Commented out `macs2`, instead loads HPC's `macs2`

- `run_snakemake.sh`: HPC has its own `snakemake` package, so no need to create one

- `run_snakemake.sh`: HPC uses a different job scheduler, so the original `bsub` command was changed to `sbatch`

- `SCENTfunctions.R`: added a `library(Matrix)` call to import package

- `pandas.yaml` and `seurat.yaml`: added a `scipy=1.14.1` or else `false_discovery_control` could not be imported

- `run_SCENT.R`: keeps giving me error messages of RuleException (caused by calculation), so I set up a check (if-else) -->