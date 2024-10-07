# source Karthik's version of SCENT (derived from SCENT GitHub)
source('workflow/scripts/SCENTfunctions.R')

library(Matrix)

# read in RNA and ATAC matrices and metadata
rna.path <-  snakemake@input[[1]]
atac.path <- snakemake@input[[2]]
metadata.path <- snakemake@input[[3]]
rna <- readRDS(rna.path)
atac <- readRDS(atac.path)
metadata <- readRDS(metadata.path)

# read in gene-peak pairs and re-order columns for SCENT
gene.peak <- read.csv(snakemake@input[[4]])
print(snakemake@input[[4]])
gene.peak <- gene.peak[ ,c('gene', 'peak')]

# create a SCENT object
scent.object <- CreateSCENTObj(
    rna = rna,
    atac = atac,
    meta.data = metadata,
    peak.info = gene.peak,
    covariates = c('log(nUMI)', 'percent.mito'),
    celltypes = 'celltype'
)

# run SCENT algorithm for user-specified cell type
scent.object <- SCENT_algorithm(
    object = scent.object,
    celltype = snakemake@params[['celltype']],
    ncores = 8
)

# write output results to a table
write.table(
    scent.object@SCENT.result,
    file = snakemake@output[[1]],
    row.names = F,
    col.names = T
)
