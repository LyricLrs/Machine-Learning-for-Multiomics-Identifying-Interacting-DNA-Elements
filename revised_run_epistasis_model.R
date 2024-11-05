# This script runs a model to test epistasis on all enhancer pairs determined
# from running SCENT.
#
# Author: Karthik Guruvayurappan

# library(Seurat) # for normalization
# library(dplyr) # for data frames
# library(boot)
library(Matrix)
# read in RNA-seq matrix, ATAC-seq matrix, and cell-level metadata
rna <- readRDS("Data/rna_matrix.rds")
atac <- readRDS("Data/atac_matrix.rds")
metadata <- readRDS("Data/metadata.rds")

# filter RNA and ATAC matrices to only include genes and peaks
# with at-least 5% non-zero counts
# rna <- rna[(rowSums(rna > 0) / ncol(rna)) >= 0.05, ]
# atac <- atac[(rowSums(atac > 0) / ncol(atac)) >= 0.05, ]

# # create Seurat object and normalize RNA-seq data
# commented out since we are now using raw count data
# print('normalizing data!')
# seurat.rna <- CreateSeuratObject(counts = rna)
# seurat.rna <- NormalizeData(seurat.rna)
# rna <- GetAssayData(seurat.rna, 'data')


# read in enhancer pairs determined from SCENT
#for testing
enhancer.pairs <- read.csv("CDS-2024-Fall-Capstone/results/enhancer_pairs/enhancer_pairs_1.csv")
enhancer.pairs <- enhancer.pairs[1,]

# create dataframe to store epistasis model results
results <- data.frame(matrix(ncol = 11, nrow = nrow(enhancer.pairs)))

# define poisson glm function (for bootstrapping)
poisson.int.coefficient <- function(data, idx = seq_len(nrow(data)), formula) {
    mdl <- glm(formula, family = 'poisson', data = data[idx,,drop = FALSE])
    mdl.values <- summary(mdl)$coefficients
    interaction.coef <- NA
    if ("enhancer.1.atac:enhancer.2.atac" %in% rownames(mdl.values)) {
        interaction.coef <- mdl.values[
            'enhancer.1.atac:enhancer.2.atac',
            'Estimate'
        ]
    }
    interaction.coef
}

## define functions
#' Interpolate a p-value from quantiles that should be "null scaled"
#'
#' @param q bootstrap quantiles, centered so that under the null, theta = 0
#' @return two-sided p-value
#' @export
interp_pval = function(q) {
  R = length(q)
  tstar = sort(q)
  zero = findInterval(0, tstar)
  if(zero == 0 || zero == R) return(2/R) # at/beyond extreme values
  pval = 2*min(zero/R, (R-zero)/R)
  pval
}


#' Derive a p-value from a vector of bootstrap samples using the "basic" calculation
#'
#' @param obs observed value of parameter (using actual data)
#' @param boot vector of bootstraps
#'
#' @return p-value
#' @export
basic_p = function(obs, boot, null = 0){
  interp_pval(2 * obs - boot - null)
}

# iterate through list of enhancer pairs and run epistasis model
for (i in 1:nrow(enhancer.pairs)) {

    # get name of enhancer 1, enhancer 2, and gene
    enhancer.1 <- enhancer.pairs$enhancer_1[i]
    enhancer.2 <- enhancer.pairs$enhancer_2[i]
    gene <- enhancer.pairs$gene[i]

    ###get the corresponding cells
    ###change the value here if need 0,1 & 1,0 & 1,1
    cells_00 <- colnames(atac)[
    atac[enhancer.1, ] == 1 & atac[enhancer.2, ] == 1
    ]

    #revised atac
    atac_00 <- atac[c(enhancer.1, enhancer.2), cells_00]

    #get all the value needs 
    enhancer.1.atac <- atac_00[enhancer.1, ]
    enhancer.2.atac <- atac_00[enhancer.2, ]
    gene.rna <- rna[gene, cells_00]
    percent.mito <- metadata[cells_00,]$percent.mito
    umis <- metadata[cells_00,]$nUMI

    #check the length of cells for each X
    print(length(gene.rna))
    print(length(enhancer.1.atac))
    print(length(enhancer.2.atac))
    print(length(percent.mito))
    print(length(umis))

    # fit linear model
    mdl <- glm(
        gene.rna ~ enhancer.1.atac * enhancer.2.atac + percent.mito + offset(log(umis)),
        family = 'poisson'
    )

    # store model values and write to data frame
    mdl.values <- summary(mdl)$coefficients

    # isolate individual values of interest
    intercept <- NA
    beta.1.estimate <- NA 
    beta.2.estimate <- NA 
    interaction.estimate <- NA 

    beta.1.pvalue <- NA 
    beta.2.pvalue <- NA 
    interaction.pvalue <- NA 

    # add error checking for possible NA values
    if ("(Intercept)" %in% rownames(mdl.values)) {
        intercept <- mdl.values['(Intercept)', 'Estimate']
    }

    if ("enhancer.1.atac" %in% rownames(mdl.values)) {
        beta.1.estimate <- mdl.values['enhancer.1.atac', 'Estimate']
        beta.1.pvalue <- mdl.values['enhancer.1.atac', 'Pr(>|z|)']
    }

    if ("enhancer.2.atac" %in% rownames(mdl.values)) {
        beta.2.estimate <- mdl.values['enhancer.2.atac', 'Estimate']
        beta.2.pvalue <- mdl.values['enhancer.2.atac', 'Pr(>|z|)']
    }

    if ("enhancer.1.atac:enhancer.2.atac" %in% rownames(mdl.values)) {
        interaction.estimate <- mdl.values[
            'enhancer.1.atac:enhancer.2.atac',
            'Estimate'
        ]
        interaction.pvalue <- mdl.values[
            'enhancer.1.atac:enhancer.2.atac',
            'Pr(>|z|)'
        ]
    }

    # implement bootstrapping
    model.df <- data.frame(cbind(gene.rna, enhancer.1.atac, enhancer.2.atac, percent.mito, umis))
    model.formula <- as.formula('gene.rna ~ enhancer.1.atac * enhancer.2.atac + percent.mito + offset(log(umis))')
    bootstrap.coefs <- boot::boot(model.df, poisson.int.coefficient, R = 100, formula = model.formula, stype = 'i', parallel = 'multicore', ncpus = 8)
    bootstrap.pvalue <- basic_p(bootstrap.coefs$t0[1], bootstrap.coefs$t[, 1])

    # implement iterative bootstrapping if p-value looks significant
    if (bootstrap.pvalue < 0.1) {
        bootstrap.coefs <- boot::boot(model.df, poisson.int.coefficient, R = 5000, formula = model.formula, stype = 'i', parallel = 'multicore', ncpus = 8)
        bootstrap.pvalue <- basic_p(bootstrap.coefs$t0[1], bootstrap.coefs$t[, 1])
    }

    # add model results to results data frame
    mdl.vector <- c(gene, enhancer.1, enhancer.2, intercept, beta.1.estimate,
                    beta.1.pvalue, beta.2.estimate, beta.2.pvalue,
                    interaction.estimate, interaction.pvalue, bootstrap.pvalue)
    results[i, ] <- mdl.vector
}

# write results to output file
write.csv(results, "model_test_output.csv", row.names = FALSE)
