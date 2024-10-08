#!/bin/bash

#SLURM -W 120:00                 # for 120 hours of wall clock time
#SLURM -J run_snakemake          # Job title
#SLURM -o run_snakemake.out   # Output file name
#SLURM -e run_snakemake.err   # Error file name

# activate snakemake mamba environment
source $HOME/.bashrc
# mamba activate snakemake

# run snakemake for whole pipeline (ending with volcano plot)
snakemake -c 4 --use-conda --jobs 32 --cluster 'SLURM -W 48:00 -n 8 -R "rusage[mem=16G]" -o out.%J.txt -e err.%J.txt' results/SCENT_peak_gene/significant_peak_gene_associations.csv

# Random comments pertaining to various components of the job submission.
# .%J adds job ID number to output files
# run using SLURM < run_snakemake.sh. '<'
