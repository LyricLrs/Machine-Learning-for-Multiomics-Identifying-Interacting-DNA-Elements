library(Matrix)

atac_matrix <- readRDS("data/PBMC/atac_matrix.rds")
rna_matrix <- readRDS("data/PBMC/rna_matrix.rds")
# metafile.rds <- readRDS("data/PBMC/metafile.rds")

atac_matrix_100 <- atac_matrix[1:100, 1:100]
rna_matrix_100 <- rna_matrix[1:100, 1:100]

saveRDS(atac_matrix_100, "data/atac_matrix.rds")
saveRDS(rna_matrix_100, "data/rna_matrix.rds")