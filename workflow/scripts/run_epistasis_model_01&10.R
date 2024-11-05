# This script runs separate models to test epistasis on enhancer pairs for cells_1 and cells_2,
# and saves the results in two separate files.
# Author: Karthik Guruvayurappan

library(Seurat) # for normalization
library(dplyr) # for data frames
library(boot)

# Read in enhancer pairs with pre-included ATAC and RNA information
# edit the paths if needed
rna <- readRDS('../../data/rna_matrix.rds')
atac <- readRDS('../../data/atac_matrix.rds')
metadata <- readRDS('../../data/meta_matrix.rds')

enhancer.pairs <- read.csv('../../data_processing/11_pairs/11_pairs2.csv')

desired_celltypes <- c('CD8-Naive')

rna <- rna[, metadata$celltype %in% desired_celltypes]
atac <- atac[, metadata$celltype %in% desired_celltypes]
metadata <- metadata[metadata$celltype %in% desired_celltypes,]

# Create dataframes to store epistasis model results for both cells_1 and cells_2
results_cells_1 <- data.frame(matrix(ncol = 6, nrow = nrow(enhancer.pairs)))
results_cells_2 <- data.frame(matrix(ncol = 6, nrow = nrow(enhancer.pairs)))

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

# Iterate over enhancer pairs and apply separate regression models for cells_1 and cells_2
for (i in 1:nrow(enhancer.pairs)) {
  enhancer.1 <- enhancer.pairs$enhancer_1[i]
  enhancer.2 <- enhancer.pairs$enhancer_2[i]
  gene <- enhancer.pairs$gene[i]
  
  # Define cell groups for cells_1 (enhancer.1 active) and cells_2 (enhancer.2 active)
  cells_1 <- colnames(atac)[atac[enhancer.2, ] == 0 & atac[enhancer.1, ] %in% c(0,1)]
  cells_2 <- colnames(atac)[atac[enhancer.1, ] == 0 & atac[enhancer.2, ] %in% c(0,1)]
  
  # Perform regression for cells_1 if any cells exist
  if (length(cells_1) > 0) {
    enhancer.atac <- atac[enhancer.1, cells_1]
    gene.rna <- rna[gene, cells_1]
    percent.mito <- metadata[cells_1,]$percent.mito
    umis <- metadata[cells_1,]$nUMI
    
    # Fit the model
    mdl <- glm(
      gene.rna ~ enhancer.atac + percent.mito + offset(log(umis)),
      family = 'poisson'
    )
    mdl.values <- summary(mdl)$coefficients
    
    # isolate individual values of interest
    intercept <- NA
    beta.estimate <- NA 
    beta.pvalue <- NA 
    
    # add error checking for possible NA values
    if ("(Intercept)" %in% rownames(mdl.values)) {
      intercept <- mdl.values['(Intercept)', 'Estimate']
    }
    
    if ("enhancer.atac" %in% rownames(mdl.values)) {
      beta.estimate <- mdl.values['enhancer.atac', 'Estimate']
      beta.pvalue <- mdl.values['enhancer.atac', 'Pr(>|z|)']
    }

    
    # Bootstrapping for cells_1
    model.df <- data.frame(cbind(gene.rna, enhancer.atac, percent.mito, umis))
    model.formula <- as.formula('gene.rna ~ enhancer.atac + percent.mito + offset(log(umis))')
    bootstrap.coefs <- boot::boot(model.df, poisson.coefficient, R = 100, formula = model.formula, stype = 'i', parallel = 'multicore', ncpus = 8)
    bootstrap.pvalue <- basic_p(bootstrap.coefs$t0[1], bootstrap.coefs$t[, 1])
    
    # implement iterative bootstrapping if p-value looks significant
    if (bootstrap.pvalue < 0.1) {
      bootstrap.coefs <- boot::boot(model.df, poisson.coefficient, R = 5000, formula = model.formula, stype = 'i', parallel = 'multicore', ncpus = 8)
      bootstrap.pvalue <- basic_p(bootstrap.coefs$t0[1], bootstrap.coefs$t[, 1])
    }
    
    # add model results to results data frame
    mdl.vector <- c(gene, enhancer.1, intercept, beta.estimate,
                    beta.pvalue, bootstrap.pvalue)
    results_cells_1[i, ] <- mdl.vector
  }
  # Perform regression for cells_2 if any cells exist
  if (length(cells_2) > 0) {
    enhancer.atac <- atac[enhancer.2, cells_2]
    gene.rna <- rna[gene, cells_2]
    percent.mito <- metadata[cells_2,]$percent.mito
    umis <- metadata[cells_2,]$nUMI
    
    # Fit the model
    mdl <- glm(
      gene.rna ~ enhancer.atac + percent.mito + offset(log(umis)),
      family = 'poisson'
    )
    mdl.values <- summary(mdl)$coefficients
    
    # isolate individual values of interest
    
    intercept <- NA
    beta.estimate <- NA 
    beta.pvalue <- NA 
    
    # add error checking for possible NA values
    if ("(Intercept)" %in% rownames(mdl.values)) {
      intercept <- mdl.values['(Intercept)', 'Estimate']
    }
    
    if ("enhancer.atac" %in% rownames(mdl.values)) {
      beta.estimate <- mdl.values['enhancer.atac', 'Estimate']
      beta.pvalue <- mdl.values['enhancer.atac', 'Pr(>|z|)']
    }
    
    # Bootstrapping for cells_2
    model.df <- data.frame(cbind(gene.rna, enhancer.atac, percent.mito, umis))
    model.formula <- as.formula('gene.rna ~ enhancer.atac + percent.mito + offset(log(umis))')
    bootstrap.coefs <- boot::boot(model.df, poisson.coefficient, R = 100, formula = model.formula, stype = 'i', parallel = 'multicore', ncpus = 8)
    bootstrap.pvalue <- basic_p(bootstrap.coefs$t0[1], bootstrap.coefs$t[, 1])
    
    # implement iterative bootstrapping if p-value looks significant
    if (bootstrap.pvalue < 0.1) {
      bootstrap.coefs <- boot::boot(model.df, poisson.coefficient, R = 5000, formula = model.formula, stype = 'i', parallel = 'multicore', ncpus = 8)
      bootstrap.pvalue <- basic_p(bootstrap.coefs$t0[1], bootstrap.coefs$t[, 1])
    }
    
    # add model results to results data frame
    mdl.vector <- c(gene, enhancer.2, intercept, beta.estimate,
                    beta.pvalue, bootstrap.pvalue)
    results_cells_2[i, ] <- mdl.vector
    
  }
}

# Write results to separate output files
write.csv(results_cells_1, '../../results/model_results_cells_1.csv', row.names = FALSE)
write.csv(results_cells_2, '../../results/model_results_cells_2.csv', row.names = FALSE)