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

3. Clone Repository to Cluster

```
git clone https://github.com/ANNIZHENG/CDS-2024-Fall-Capstone.git
```

4. Request Space (1 hour, 4 core, 100G memory)

```
srun -t 1:00:00 -c 4 --mem=100G --pty /bin/bash
```

5. Run Snakemake in HPC
```
bash run_snakemake.sh
```

6. Some Useful Commands
```
# Remove all created snakemake files (clear memory)
rm -rf .snakemake/conda/*

# Remove the cache file (clear memory)
rm -rf cache

# Remove any produced file (clear memory)
rm  -rf *.txt
rm -rf results

# Update the console's GitHub reposotiry
git pull origin main
```
<!-- rm -rf resources -->

## Modifications

- `seurat.yaml`: Commented out `macs2`, instead loads HPC's `macs2`

- `run_snakemake.sh`: HPC uses a different job scheduler, so the original `bsub` command was changed to `sbatch`

- `run_snakemake.sh`: HPC has its own `snakemake` package <-- we would use this. 

- `SCENTfunctions.R`: added a `library(Matrix)` call to import package

- `seurat.yaml`: added a `scipy=1.11.1` or else `false_discovery_control` could not be imported

- `run_SCENT.R` and ``: keeps giving me error messages of RuleException (caused by calculation), so I set up a check (if-else)

Below is the error message. 
I temporarily prevent it with an if-else check on `nonzero...` variables. 

```
Warning message:
undefined slot classes in definition of "SCENT": rna(class "dgCMatrix"), atac(class "dgCMatrix") 
Error in if (nonzero_m > 0.05 & nonzero_a > 0.05) { : 
  missing value where TRUE/FALSE needed
Calls: SCENT_algorithm

RuleException:
CalledProcessError in line 83 of /scratch/.../workflow/Snakefile:
Command 'source /share/apps/mambaforge/23.1.0/bin/activate '/scratch/.../.snakemake/conda/b86
91cb1389694139f08589ae49b38cc'; set -euo pipefail;  Rscript --vanilla /scratch/.../.snakemake
/scripts/tmp3s56enug.run_SCENT.R' returned non-zero exit status 1.
  File "/scratch/az1932/CDS-2024-Fall-Capstone/workflow/Snakefile", line 83, in __rule_run_SCENT
  File "/share/apps/python/3.8.6/intel/lib/python3.8/concurrent/futures/thread.py", line 57, in run
```