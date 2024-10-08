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

2. Open Cluster

3. Clone Repository to Cluster
```console
git clone https://github.com/ANNIZHENG/CDS-2024-Fall-Capstone.git
```

4. Request Space
```
srun -t 1:00:00 -c 4 --mem=16000 --pty /bin/bash
```

5. Snakemake Set Up

```console
cd CDS-2024-Fall-Capstone
mkdir -p cache
export CONDA_PKGS_DIRS=cache
module load mambaforge/23.1.0
module load snakemake/6.12.3
module load macs2/intel/2.2.7.1
```
<!-- snakemake -c 4 -->

6. Run Model 
```
bash run_snakemake.sh
```