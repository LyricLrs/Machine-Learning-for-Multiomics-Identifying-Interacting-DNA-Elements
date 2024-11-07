# This script runs a single model to test epistasis on enhancer pairs for combined cells_1 and cells_2.
# Author: Karthik Guruvayurappan

library(Seurat) # for normalization
library(dplyr) # for data frames
library(boot)

# Read in enhancer pairs with pre-included ATAC and RNA information
# edit the paths if needed
rna <- readRDS(snakemake@input[[1]])
atac <- readRDS(snakemake@input[[2]])
metadata <- readRDS(snakemake@input[[3]])
desired_celltypes <- snakemake@params[['celltype']]

# Combine all data_processing/11_pairs{i}.csv files
enhancer.pairs <- data.frame()
for (file in snakemake@input[4:length(snakemake@input)]) {
    enhancer_temp.pairs <- read.csv(file)
    enhancer.pairs <- rbind(enhancer.pairs, enhancer_temp.pairs)
}
rna <- rna[, metadata$celltype %in% desired_celltypes]
atac <- atac[, metadata$celltype %in% desired_celltypes]
metadata <- metadata[metadata$celltype %in% desired_celltypes,]

# Create a dataframe to store combined model results
results <- data.frame(matrix(ncol = 6, nrow = nrow(enhancer.pairs)))
colnames(results) <- c("gene", "enhancer", "intercept", "beta.estimate", "beta.pvalue", "bootstrap.pvalue")

# Define Poisson GLM function for bootstrapping
poisson.coefficient <- function(data, idx = seq_len(nrow(data)), formula) {
  mdl <- glm(formula, family = 'poisson', data = data[idx,,drop = FALSE])
  mdl.values <- summary(mdl)$coefficients
  interaction.coef <- NA
  if ("enhancer.atac" %in% rownames(mdl.values)) {
    interaction.coef <- mdl.values['enhancer.atac', 'Estimate']
  }
  interaction.coef
}

# Define p-value calculation functions
interp_pval <- function(q) {
  R <- length(q)
  tstar <- sort(q)
  zero <- findInterval(0, tstar)
  if(zero == 0 || zero == R) return(2/R) 
  pval <- 2*min(zero/R, (R-zero)/R)
  pval
}

basic_p <- function(obs, boot, null = 0) {
  interp_pval(2 * obs - boot - null)
}

# Iterate over enhancer pairs and apply a single combined regression model for cells_1 and cells_2
for (i in 1:nrow(enhancer.pairs)) {
  enhancer.1 <- enhancer.pairs$enhancer_1[i]
  enhancer.2 <- enhancer.pairs$enhancer_2[i]
  gene <- enhancer.pairs$gene[i]
  
  # Define cell groups for cells_1 (enhancer.1 active) and cells_2 (enhancer.2 active)
  cells_1 <- colnames(atac)[atac[enhancer.1, ] == 0 & atac[enhancer.2, ] %in% c(0,1)]
  cells_2 <- colnames(atac)[atac[enhancer.2, ] == 0 & atac[enhancer.1, ] %in% c(0,1)]
  
  # Combine cells_1 and cells_2 into one dataset
  combined_cells <- unique(c(cells_1, cells_2))
  
  if (length(combined_cells) > 0) {
    # Use combined cells for regression
    enhancer.atac <- atac[enhancer.1, combined_cells] + atac[enhancer.2, combined_cells]
    gene.rna <- rna[gene, combined_cells]
    percent.mito <- metadata[combined_cells,]$percent.mito
    umis <- metadata[combined_cells,]$nUMI
    
    # Fit the model with combined data
    mdl <- glm(
      gene.rna ~ enhancer.atac + percent.mito + offset(log(umis)),
      family = 'poisson'
    )
    mdl.values <- summary(mdl)$coefficients
    
    # Isolate individual values of interest
    intercept <- if ("(Intercept)" %in% rownames(mdl.values)) mdl.values['(Intercept)', 'Estimate'] else NA
    beta.estimate <- if ("enhancer.atac" %in% rownames(mdl.values)) mdl.values['enhancer.atac', 'Estimate'] else NA
    beta.pvalue <- if ("enhancer.atac" %in% rownames(mdl.values)) mdl.values['enhancer.atac', 'Pr(>|z|)'] else NA
    
    # Bootstrapping for combined cells
    model.df <- data.frame(cbind(gene.rna, enhancer.atac, percent.mito, umis))
    model.formula <- as.formula('gene.rna ~ enhancer.atac + percent.mito + offset(log(umis))')
    bootstrap.coefs <- boot::boot(model.df, poisson.coefficient, R = 100, formula = model.formula, stype = 'i', parallel = 'multicore', ncpus = 8)
    bootstrap.pvalue <- basic_p(bootstrap.coefs$t0[1], bootstrap.coefs$t[, 1])
    
    # Implement iterative bootstrapping if p-value looks significant
    if (bootstrap.pvalue < 0.1) {
      bootstrap.coefs <- boot::boot(model.df, poisson.coefficient, R = 5000, formula = model.formula, stype = 'i', parallel = 'multicore', ncpus = 8)
      bootstrap.pvalue <- basic_p(bootstrap.coefs$t0[1], bootstrap.coefs$t[, 1])
    }
    
    # Add model results to results data frame
    mdl.vector <- c(gene, paste(enhancer.1, enhancer.2, sep = "_"), intercept, beta.estimate, beta.pvalue, bootstrap.pvalue)
    results[i, ] <- mdl.vector
  }
}

# Write results to a single output file
write.csv(results_cells_1, snakemake@output[[1]], row.names = FALSE)