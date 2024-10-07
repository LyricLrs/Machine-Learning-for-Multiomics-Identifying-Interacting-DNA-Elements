# This script reads in expression and peak matrices from a 10X multiome
# experiment, isolates the gene and peak names, and writes them to output CSV
# files.
#
# input[[1]]: file path to expression matrix (raw counts)
# input[[2]]: file path to peak matrix (binarized to 0 to 1)
# Author: Karthik Guruvayurappan

# load 'Matrix' library
library(Matrix)

# get filenames for RNA and ATAC matrices
rna.path <- snakemake@input[[1]]
atac.path <- snakemake@input[[2]]

# read in RDS structures for RNA and ATAC matrices
rna <- readRDS(rna.path)
atac <- readRDS(atac.path)

# filter rna and atac matrices to have minimum 5% non-zero counts for each gene
rna <- rna[(rowSums(rna > 0) / ncol(rna)) >= 0.05, ]
atac <- atac[(rowSums(atac > 0) / ncol(atac)) >= 0.05, ]

# save gene and peak names to vectors
genes <- rownames(rna)
peaks <- rownames(atac)

# write gene and peak names to output CSV files
write.csv(genes, snakemake@output[[1]], row.names = FALSE)
write.csv(peaks, snakemake@output[[2]], row.names = FALSE)
