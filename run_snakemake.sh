#!/bin/bash

#SBATCH -J run_snakemake          # Job title
#SBATCH -o run_snakemake.out      # Output file name
#SBATCH -e run_snakemake.err      # Error file name

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
# snakemake --forceall --use-conda --jobs 32 --cluster 'sbatch --time=48:00:00 --mem=16G --cpus-per-task=8 -o out.%J.txt -e err.%J.txt'
snakemake --forceall --use-conda --jobs 32 --cluster 'sbatch --time=48:00:00 --mem=16G --cpus-per-task=8 -o out.%J.txt -e err.%J.txt' -p \
$(echo results/epistasis_models/epistasis_models_01_10_1_{1..32}.csv results/epistasis_models/epistasis_models_01_10_2_{1..32}.csv results/epistasis_models/epistasis_models_11_00_{1..32}.csv)

# Random comments pertaining to various components of the job submission.
# .%J adds job ID number to output files
# run using bsub < run_snakemake.sh. '<'
