# This script generates 3 different tables of significant epistatic models
# using 3 different significance thresholds. The first significance threshold
# is an FDR < 0.1 across all possible pairs tested (most stringent). The
# second threshold is a gene-wise threshold of FDR < 0.01. The third threshold
# does not use FDR correction, and simply takes all of the models with a
# analytical p-value below 5e-4.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

# declare number of batches as constant
NUM_BATCHES = 32

# concatenate results across all epistatic models
epistasis_models = pd.DataFrame()

for i in range(NUM_BATCHES):
    filename = snakemake.input[i]
    batch_models = pd.read_csv(filename)
    epistasis_models = pd.concat([epistasis_models, batch_models])

epistasis_models.columns = [
    'gene', 'enhancer_1', 'enhancer_2', 'intercept', 'beta_1_estimate',
    'beta_1_pvalue', 'beta_2_estimate', 'beta_2_pvalue',
    'interaction_estimate', 'interaction_pvalue', 'bootstrap_pvalue'
]

# compute significant results for threshold 3
significant_results = epistasis_models[
    epistasis_models['bootstrap_pvalue'] < 0.05
]

# write output to CSV file
significant_results.to_csv(
    snakemake.output[0],
    index=False
)

# compute significant results for threshold 1
epistasis_models = epistasis_models.dropna()
epistasis_models['fdr_pvalue'] = false_discovery_control(
    epistasis_models['bootstrap_pvalue']
)

significant_results = epistasis_models[
    epistasis_models['fdr_pvalue'] < 0.1
]

# write output to CSV file
significant_results.to_csv(
    snakemake.output[1],
    index=False
)

# filter for significant results from threshold 2
gene_adj_epistasis_models = pd.DataFrame()

genes = epistasis_models['gene'].unique()

for gene in genes:

    # get models for that gene
    gene_models = epistasis_models.groupby('gene').get_group(gene)

    # perform FDR correction on p-values for gene
    gene_models['fdr_pvalue'] = false_discovery_control(
        gene_models['bootstrap_pvalue']
    )

    # add gene-adjusted p-values to combined data frame
    gene_adj_epistasis_models = (
        pd.concat([gene_adj_epistasis_models, gene_models])
    )

# write significant results to dataframe
epistasis_models = gene_adj_epistasis_models
significant_results = epistasis_models[
    epistasis_models['fdr_pvalue'] < 0.1
]

# write output to CSV file
significant_results.to_csv(
    snakemake.output[2],
    index=False
)
