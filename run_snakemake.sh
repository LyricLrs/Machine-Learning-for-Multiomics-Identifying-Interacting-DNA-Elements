#!/bin/bash

#BSUB -W 120:00                 # for 120 hours of wall clock time
#BSUB -J run_snakemake          # Job title
#BSUB -o run_snakemake.out      # Output file name
#BSUB -e run_snakemake.err      # Error file name

# Load necessary modules
module load mambaforge/23.1.0
module load snakemake/6.12.3
module load macs2/intel/2.2.7.1

# create local cache
mkdir -p conda_cache/envs
mkdir -p conda_cache/pkgs
export CONDA_ENVS_PATH=conda_cache/envs
export CONDA_PKGS_DIRS=conda_cache/pkgs

# Unlock Snakemake in case of a stale lock
snakemake --unlock

# Run Snakemake with Conda, specifying your local Conda cache and cluster submission settings
snakemake --use-conda --conda-prefix conda_cache/envs \
    --jobs 32 \
    --cluster 'bsub -W 48:00 -n 8 -R "rusage[mem=16G]" -o out.%J.txt -e err.%J.txt'