# This script filters gene-peak pairs output by bedtools to only include genes
# with RNA-seq data.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# read in gene-peak pairs from bedtools
bed_names = [
    'gene_chrom',
    'gene_start',
    'gene_end',
    'gene',
    'peak_chrom',
    'peak_start',
    'peak_end'
]
gene_peak_pairs = pd.read_csv(
    snakemake.input[0],
    sep='\t',
    names=bed_names

)

# read in gene names from RNA-seq matrix
genes = pd.read_csv(snakemake.input[1])

# filter gene-peak pairs
gene_peak_pairs = gene_peak_pairs[gene_peak_pairs['gene'].isin(genes['x'])]

# adjust back to GTF format for peak starts
gene_peak_pairs['peak_start'] = gene_peak_pairs['peak_start'] + 1

# format file with only genes and peaks (for SCENT input)
gene_peak_pairs['peak_start'] = gene_peak_pairs['peak_start'].astype(str)
gene_peak_pairs['peak_end'] = gene_peak_pairs['peak_end'].astype(str)
gene_peak_pairs['peak'] = (gene_peak_pairs['peak_chrom'] + 
                           '-' +
                           gene_peak_pairs['peak_start'] +
                           '-' +
                           gene_peak_pairs['peak_end'])
gene_peak_pairs = gene_peak_pairs[['peak', 'gene']]

# read in gene names for highly variable features
genes = pd.read_csv(snakemake.input[2])

# filter for only gene-peak pairs in highly variable features
gene_peak_pairs = gene_peak_pairs[gene_peak_pairs['gene'].isin(genes['x'])]

# split gene-peak pairs into 32 dataframes (parallelization)
gene_peak_split_pairs = np.array_split(gene_peak_pairs, 32)

# write files to output CSV files
for i in np.arange(32):
    gene_peak_split_pairs[i].to_csv(
        snakemake.output[i],
        index=False)
