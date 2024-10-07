"""Annotate Gene and Peak Names

This script annotates a set of genes and peaks with name, start, and end
information and outputs .bed files for compatiblity with bedtools.

Author: Karthik Guruvayurappan
"""

# import numpy and pandas (working with DataFrames)
import numpy as np
import pandas as pd

# read in Gencode annotations
gencode_names = [
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'
    ]
annotations = pd.read_csv(
    snakemake.input[2],
    skiprows=5,
    sep='\t',
    names=gencode_names)

# filter to only include genes
annotations = annotations[annotations['feature'] == 'gene']
# get gene names from the attribute column of the GTF file
def get_gene_name(attribute):
    attribute_info = attribute.split(';')
    return attribute_info[2].split('"')[1]

gene_names = annotations['attribute'].apply(get_gene_name)
annotations['gene_name'] = gene_names

# get gene name, chromosome, start, and end
annotations = annotations[['seqname', 'start', 'end', 'gene_name']]
annotations.columns = ['chrom', 'start', 'end', 'gene']

# adjust GTF file for .bed format and save as a .bed file
annotations['start'] = annotations['start'] - 1
annotations.to_csv(
    snakemake.output[0],
    sep='\t',
    index=False,
    header=False)

# read in peaks file
peaks = pd.read_csv(snakemake.input[1])
peaks.columns = ['peak']

# parse peaks into a 3 column data frame
peaks_list = list(peaks['peak'].str.split('-'))
peaks = pd.DataFrame(peaks_list, columns=['chrom', 'start', 'end'])

# adjust dataframe to match .bed format from GTF format
peaks['start'] = peaks['start'].astype(np.int64)
peaks['end'] = peaks['end'].astype(np.int64)
peaks['start'] = peaks['start'] - 1

# write out ATAC-seq peaks as a .bed file
peaks.to_csv(snakemake.output[1], sep='\t', index=False, header=False)
