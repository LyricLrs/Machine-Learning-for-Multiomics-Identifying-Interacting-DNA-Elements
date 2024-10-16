library(Matrix)  

atac <- readRDS("Data/atac.rds")
subset_data <- atac[1:100, ]  # Adjust this range as needed
print(subset_data)

cell_names <- colnames(subset_data)
peak_pairs <- combn(nrow(subset_data), 2)

zero_zero_list <- list()  # Peak1 == 0 & Peak2 == 0
zero_one_list  <- list()  # Peak1 == 0 & Peak2 == 1
one_zero_list  <- list()  # Peak1 == 1 & Peak2 == 0
one_one_list   <- list()  # Peak1 == 1 & Peak2 == 1

for (i in 1:ncol(peak_pairs)) {
  peak1 <- as.matrix(subset_data[peak_pairs[1, i], ])
  peak2 <- as.matrix(subset_data[peak_pairs[2, i], ])
  
  # Find the cells corresponding to each combination:
  zero_zero_cells <- which(peak1 == 0 & peak2 == 0)  # peak1 == 0 and peak2 == 0
  zero_one_cells  <- which(peak1 == 0 & peak2 == 1)  # peak1 == 0 and peak2 == 1
  one_zero_cells  <- which(peak1 == 1 & peak2 == 0)  # peak1 == 1 and peak2 == 0
  one_one_cells   <- which(peak1 == 1 & peak2 == 1)  # peak1 == 1 and peak2 == 1

  # Store the values and corresponding cell names as a list for each condition
  if (length(zero_zero_cells) > 0) {
    zero_zero_list[[i]] <- list(Values = subset_data[c(peak_pairs[1, i], peak_pairs[2, i]), zero_zero_cells, drop=FALSE],
                                CellNames = cell_names[zero_zero_cells])
  }
  if (length(zero_one_cells) > 0) {
    zero_one_list[[i]] <- list(Values = subset_data[c(peak_pairs[1, i], peak_pairs[2, i]), zero_one_cells, drop=FALSE],
                               CellNames = cell_names[zero_one_cells])
  }
  if (length(one_zero_cells) > 0) {
    one_zero_list[[i]] <- list(Values = subset_data[c(peak_pairs[1, i], peak_pairs[2, i]), one_zero_cells, drop=FALSE],
                               CellNames = cell_names[one_zero_cells])
  }
  if (length(one_one_cells) > 0) {
    one_one_list[[i]] <- list(Values = subset_data[c(peak_pairs[1, i], peak_pairs[2, i]), one_one_cells, drop=FALSE],
                              CellNames = cell_names[one_one_cells])
  }
}

zero_zero_list[[1]]

#store data
saveRDS(zero_zero_list, "zero_zero_list_with_cells.rds")
saveRDS(zero_one_list, "zero_one_list_with_cells.rds")
saveRDS(one_zero_list, "one_zero_list_with_cells.rds")
saveRDS(one_one_list, "one_one_list_with_cells.rds")



#How do we find the value for a particular gene
gene_of_interest <- "AAAGCCCGTTATCCTA-1"
gene_index <- which(cell_names == gene_of_interest)

pair_1_values <- zero_zero_list[[1]]$Values  # Get the values matrix
pair_1_gene_value <- pair_1_values[, gene_index, drop=FALSE]  # Extract the values for the gene

pair_1_gene_value
