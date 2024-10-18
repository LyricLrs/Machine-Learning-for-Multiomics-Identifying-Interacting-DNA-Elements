# install.packages("Matrix")
# install.packages('IRkernel')
library(Matrix)

atac <- readRDS("data/atac_matrix.rds")

#will be revised later
subset_data <- atac[1:100, 1:10] 
print(subset_data)

peak_names <- rownames(subset_data)
gene_names <- colnames(subset_data)

gene_zero_zero_list <- list()
gene_zero_one_list <- list()
gene_one_zero_list <- list()
gene_one_one_list <- list()

output_list <- list()
for (gene_index in 1:ncol(subset_data)) {
  gene_values <- subset_data[, gene_index]
  
  peak_pairs <- combn(nrow(subset_data), 2)
  
  zero_zero_pairs <- list()
  zero_one_pairs <- list()
  one_zero_pairs <- list()
  one_one_pairs <- list()

  for (i in 1:ncol(peak_pairs)) {
    peak1_value <- gene_values[peak_pairs[1, i]]
    peak2_value <- gene_values[peak_pairs[2, i]]

    peak1_name <- peak_names[peak_pairs[1, i]]  
    peak2_name <- peak_names[peak_pairs[2, i]]  

    #0 0 
    if (peak1_value == 0 & peak2_value == 0) {
      zero_zero_pairs[[length(zero_zero_pairs) + 1]] <- paste(gene_names[gene_index], peak1_name, peak2_name, sep = ",")
    }
    
    #0 1
    if (peak1_value == 0 & peak2_value == 1) {
      zero_one_pairs[[length(zero_one_pairs) + 1]] <- paste(gene_names[gene_index], peak1_name, peak2_name, sep = ",")
    }
  
    #1 0 
    if (peak1_value == 1 & peak2_value == 0) {
      one_zero_pairs[[length(one_zero_pairs) + 1]] <- paste(gene_names[gene_index], peak1_name, peak2_name, sep = ",")
    }
    
    #1 1
    if (peak1_value == 1 & peak2_value == 1) {
      one_one_pairs[[length(one_one_pairs) + 1]] <- paste(gene_names[gene_index], peak1_name, peak2_name, sep = ",")
    }

    }
  gene_zero_zero_list[[gene_names[gene_index]]] <- zero_zero_pairs
  gene_zero_one_list[[gene_names[gene_index]]] <- zero_one_pairs
  gene_one_zero_list[[gene_names[gene_index]]] <- one_zero_pairs
  gene_one_one_list[[gene_names[gene_index]]] <- one_one_pairs
  }


gene_zero_zero_list <- unlist(gene_zero_zero_list)
gene_zero_one_list <- unlist(gene_zero_one_list)
gene_one_zero_list <- unlist(gene_one_zero_list)
gene_one_one_list <- unlist(gene_one_one_list)

if (!dir.exists("Processed_data_example")) {
  dir.create("Processed_data_example")
}

writeLines(gene_zero_zero_list, "Processed_data_example/gene_zero_zero_list.csv")
writeLines(gene_zero_one_list, "Processed_data_example/gene_zero_one_list.csv")
writeLines(gene_one_zero_list, "Processed_data_example/gene_one_zero_list.csv")
writeLines(gene_one_one_list, "Processed_data_example/gene_one_one_list.csv")




