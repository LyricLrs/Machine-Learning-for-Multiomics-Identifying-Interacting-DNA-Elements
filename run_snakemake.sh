#!/bin/bash

#BSUB -W 120:00                 # for 120 hours of wall clock time
#BSUB -J run_snakemake          # Job title
#BSUB -o run_snakemake.out   # Output file name
#BSUB -e run_snakemake.err   # Error file name

# To accommodate HPC's environment
mkdir -p cache
export CONDA_PKGS_DIRS=cache
module load mambaforge/23.1.0
module load snakemake/6.12.3
module load macs2/intel/2.2.7.1

# activate snakemake mamba environment
source $HOME/.bashrc

# HPC has its own snakemake so we would use that
# mamba activate snakemake

# run snakemake for whole pipeline (ending with volcano plot)
# snakemake --use-conda --jobs 32 --cluster 'bsub -W 48:00 -n 8 -R "rusage[mem=16G]" -o out.%J.txt -e err.%J.txt' results/SCENT_peak_gene/significant_peak_gene_associations.csv
snakemake -c 4 --use-conda --jobs 32 --cluster 'sbatch --time=48:00:00 --cpus-per-task=8 --mem=16G -o out.%J.txt -e err.%J.txt' results/SCENT_peak_gene/significant_peak_gene_associations.csv


# Random comments pertaining to various components of the job submission.
# .%J adds job ID number to output files
# run using bsub < run_snakemake.sh. '<'
