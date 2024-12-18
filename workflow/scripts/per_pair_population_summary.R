library(Seurat) 
library(dplyr) # for data frames
library(boot)
library(ggplot2)

library(pscl)
# edit the paths if needed
rna <- readRDS('/Users/huajingru/Desktop/Fall_2024/Capstone/PBMC/rna_matrix.rds')
atac <- readRDS('/Users/huajingru/Desktop/Fall_2024/Capstone/PBMC/atac_matrix.rds')
metadata <- readRDS('/Users/huajingru/Desktop/Fall_2024/Capstone/PBMC/metafile.rds')
#enhancer.pairs <- read.csv('/Users/huajingru/Desktop/Fall_2024/Capstone/CD14-Mono/enhancer_pairs_CD14-Mono/enhancer_pairs_1.csv')

directory_path <- "/Users/huajingru/Desktop/Fall_2024/Capstone/CD14-Mono/enhancer_pairs_CD14-Mono/"
csv_files <- list.files(directory_path, pattern = "\\.csv$", full.names= TRUE)
enhancer_data_list <- lapply(csv_files, read.csv)
enhancer.pairs <- do.call(rbind, enhancer_data_list)

desired_celltypes <- c('CD14-Mono')

rna <- rna[, metadata$celltype %in% desired_celltypes]
atac <- atac[, metadata$celltype %in% desired_celltypes]
metadata <- metadata[metadata$celltype %in% desired_celltypes,]

# Initialize a data frame to store per-pair population statistics
per_pair_summary <- data.frame(
  enhancer_1 = character(),
  enhancer_2 = character(),
  gene = character(),
  population = character(),
  total_cells = numeric(),
  zero_expression_count = numeric(),
  zero_expression_percentages = numeric(),
  stringsAsFactors = FALSE
)

# Loop through enhancer pairs
for (i in 1:nrow(enhancer.pairs)) {

  enhancer.1 <- enhancer.pairs$enhancer_1[i]
  enhancer.2 <- enhancer.pairs$enhancer_2[i]
  gene <- enhancer.pairs$gene[i]
  
  # Define cells in each group based on conditions
  cells_00 <- colnames(atac)[atac[enhancer.1, ] == 0 & atac[enhancer.2, ] == 0]
  cells_01 <- colnames(atac)[atac[enhancer.1, ] == 0 & atac[enhancer.2, ] == 1]
  cells_10 <- colnames(atac)[atac[enhancer.1, ] == 1 & atac[enhancer.2, ] == 0]
  cells_11 <- colnames(atac)[atac[enhancer.1, ] == 1 & atac[enhancer.2, ] == 1]
  
  # Define expression data for the gene in each group
  gene_00 <- rna[gene, cells_00]
  gene_01 <- rna[gene, cells_01]
  gene_10 <- rna[gene, cells_10]
  gene_11 <- rna[gene, cells_11]
  
  # Create a temporary data frame for this enhancer-gene pair
  pair_summary <- data.frame(
    enhancer_1 = enhancer.1,
    enhancer_2 = enhancer.2,
    gene = gene,
    population = c("00", "01", "10", "11"),
    total_cells = c(length(cells_00), length(cells_01), length(cells_10), length(cells_11)),
    zero_expression_count = c(sum(gene_00 == 0), sum(gene_01 == 0), sum(gene_10 == 0), sum(gene_11 == 0)),
    stringsAsFactors = FALSE
  )
  
  # Calculate zero expression fractions
  pair_summary$zero_expression_percentages <- pair_summary$zero_expression_count / pair_summary$total_cells
  
  # Append to the main per-pair summary
  per_pair_summary <- rbind(per_pair_summary, pair_summary)
}

# Print and optionally save the per-pair summary
cat("Per-Pair Population Summary:\n")
print(head(per_pair_summary))
write.csv(per_pair_summary, "/Users/huajingru/Desktop/Fall_2024/Capstone/CD14-Mono/per_pair_population_summary.csv", row.names = FALSE)
