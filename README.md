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

```
# Create
git clone https://github.com/ANNIZHENG/CDS-2024-Fall-Capstone.git

# Update
git pull origin main
```

4. Send task to HPC Cluster
```
sbatch run_snakemake.sh
```

5. Useful commands to clear some memory space

```
rm -rf *.txt
rm -rf *.out
rm -rf cache
rm -rf .snakemake/conda/
conda clean --all
```

## Modifications

- `seurat.yaml`: Commented out `macs2`, instead loads HPC's `macs2`

- `run_snakemake.sh`: HPC has its own `snakemake` package, so no need to create one

- `run_snakemake.sh`: HPC uses a different job scheduler, so the original `bsub` command was changed to `sbatch`

- `SCENTfunctions.R`: added a `library(Matrix)` call to import package

- `pandas.yaml` and `seurat.yaml`: added a `scipy=1.14.1` or else `false_discovery_control` could not be imported

**Import Error**
```
Traceback (most recent call last):
  File ".snakemake/scripts/tmpv9ttqe5u.find_enhancer_pairs.py", line 13, in <module>
    from scipy.stats import false_discovery_control
  
ImportError: cannot import name 'false_discovery_control' from 'scipy.stats' (/share/apps/python/3.8.6/intel/lib/python3.8/site-packages/scipy-1.5.2-py3.8-linux-x86_64.egg/scipy/stats/__init__.py)
```

- `run_SCENT.R`: keeps giving me error messages of RuleException (caused by calculation), so I set up a check (if-else)

**Rule Error**
```
Warning message:
undefined slot classes in definition of "SCENT": rna(class "dgCMatrix"), atac(class "dgCMatrix") 
Error in if (nonzero_m > 0.05 & nonzero_a > 0.05) { : 
  missing value where TRUE/FALSE needed
Calls: SCENT_algorithm
```