# This scripts concatenates separate output files from SCENT and determines all
# possible enhancer pairs that can be used for testing epistasis.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd
import itertools
## from scipy.stats import false_discovery_control
from statsmodels.stats.multitest import multipletests


# concatenate all peak-gene results across 16 SCENT batch runs
peak_gene_results = []

for i in np.arange(32):
    peak_gene_results.append(
        pd.read_csv(snakemake.input[i], sep=' ')
    )
peak_gene_results = pd.concat(peak_gene_results)

# for i in np.arange(32):
#     peak_gene_results.append(
#         pd.read_csv('/data/deyk/karthik/multiome_epistasis/results/SCENT_outputs/SCENT_output_' + str(i+1) + '.csv', sep=' ')
#     )
# peak_gene_results = pd.concat(peak_gene_results)

# perform multiple testing correction on p-values
## peak_gene_results['adj_p'] = false_discovery_control(peak_gene_results['boot_basic_p'])
reject, pvals_corrected, _, _ = multipletests(peak_gene_results['boot_basic_p'], alpha=0.05, method='fdr_bh')
peak_gene_results['adj_p'] = pvals_corrected

# filter for all enhancer-gene pairs with an FDR < 0.1
peak_gene_results = peak_gene_results[peak_gene_results['adj_p'] < 0.1]

# write significant peak-gene pairs to output file (for reference)
peak_gene_results.to_csv(
    snakemake.output[0],
    index = False
)

# write significant genes to separate file (for GSEA)
unique_genes = peak_gene_results['gene'].unique()
pd.Series(unique_genes).to_csv(snakemake.output[33], index = False)

# get list of genes that have more than 1 significant peak-gene link
gene_peak_counts = peak_gene_results['gene'].value_counts()
multi_peak_genes = pd.Series(gene_peak_counts.index[gene_peak_counts > 1])

# filter for peak-gene pairs where genes have multiple linked peaks
multi_peak_peaks = peak_gene_results['gene'].isin(multi_peak_genes)
peak_gene_results = peak_gene_results[multi_peak_peaks]

# determine all possible enhancer pairs for the same gene
enhancer_pair_lists = []

for gene in multi_peak_genes:

    # determine pairs of peaks
    gene_peaks = peak_gene_results[peak_gene_results['gene'] == gene]
    peaks = gene_peaks['peak']
    enhancer_pairs = itertools.combinations(peaks, 2)
    enhancer_pair_lists.append(list(enhancer_pairs))

# create dataframe with all enhancer pairs
enhancer_pair_df = pd.DataFrame()
enhancer_pair_df['gene'] = multi_peak_genes
enhancer_pair_df['enhancer_pairs'] = enhancer_pair_lists
enhancer_pair_df = enhancer_pair_df.explode('enhancer_pairs')

# separate lists into separate dataframe columns for each enhancer
enhancer_1 =  enhancer_pair_df['enhancer_pairs'].apply(lambda x: x[0])
enhancer_2 = enhancer_pair_df['enhancer_pairs'].apply(lambda x: x[1])
enhancer_pair_df['enhancer_1'] = enhancer_1
enhancer_pair_df['enhancer_2'] = enhancer_2
enhancer_pair_df = enhancer_pair_df[['gene', 'enhancer_1', 'enhancer_2']]

# split enhancer pairs into 32 dataframes (parallelization)
enhancer_pair_split_df = np.array_split(enhancer_pair_df, 32)

# write files to output CSV files
for i in np.arange(32):
    enhancer_pair_split_df[i].to_csv(
        snakemake.output[i+1],
        index=False)
