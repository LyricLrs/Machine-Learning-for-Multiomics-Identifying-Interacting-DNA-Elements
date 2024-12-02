# This script identifies highly variable genes for epistasis testing using the
# FindVariableFeatures functionality (mvp method) implemented in Seurat. This
# method first computes the mean and dispersion using the raw counts for each
# gene. Then, by binning the genes into 20 bins, this method identifies
# outliers for dispersion based on z-score, and uses those as highly variable
# genes for downstream single-cell analysis.
#
# Author: Karthik Guruvayurappan

library(Seurat)

# read in scRNA-seq matrix into Seurat
rna <- readRDS(snakemake@input[[1]])
metadata <- readRDS(snakemake@input[[2]])

desired_celltypes <- snakemake@params[['celltype']]

rna <- rna[, metadata$celltype %in% desired_celltypes]

# rna <- CreateSeuratObject(counts = rna)

# compute highly variable features using the SCTransform method
# rna <- SCTransform(rna)

# save highly variable genes to an output CSV file
# var.genes <- rna@assays$SCT@var.features
var.genes <- rownames(rna)
write.csv(
    var.genes,
    snakemake@output[[1]],
    row.names = FALSE,
    quote = FALSE
)
