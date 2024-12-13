# This script plots an interaction term which summarizes the interaction terms
# observed from our epistasis models.
#
# Author: Karthik Guruvayurappan

library(stats)
library(ggplot2) # for plotting
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)

# iterate through different model output files
epistasis.models <- data.frame()

#enhancer1
folder_path <- "/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/results/epistasis_model_output/CD14-mono/enhancer1/"
file_list <- list.files(path = folder_path, pattern = "epistasis_models_01_10_1_.*\\.csv", full.names = TRUE)
read_custom_csv <- function(file) {
  read.csv(file, col.names = c("gene", "enhancer1", "enhancer2", "intercept", "beta_estimate", "beta_pvalue", "bootstrap_pvalue"))
}
enhancer1 <- do.call(rbind, lapply(file_list, read_custom_csv))
head(enhancer1)

#enhancer2
folder_path <- "/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/results/epistasis_model_output/CD14-mono/enhancer2/"
file_list <- list.files(path = folder_path, pattern = "epistasis_models_01_10_2_.*\\.csv", full.names = TRUE)
read_custom_csv <- function(file) {
  read.csv(file, col.names = c("gene", "enhancer1", "enhancer2", "intercept", "beta_estimate", "beta_pvalue", "bootstrap_pvalue"))
}
enhancer2 <- do.call(rbind, lapply(file_list, read_custom_csv))
head(enhancer2)

#interaction
folder_path <- "/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/results/epistasis_model_output/CD14-mono/interaction/"
file_list <- list.files(path = folder_path, pattern = "epistasis_models_11_00_.*\\.csv", full.names = TRUE)
read_custom_csv <- function(file) {
  read.csv(file, col.names = c("gene", "enhancers", "intercept", "beta.estimate", "beta.pvalue", "bootstrap.pvalue"))
}
interaction <- do.call(rbind, lapply(file_list, read_custom_csv))
head(interaction)



# # read in epistasis models from batch
# enhancer1 <- read.csv('/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/results/hpc_model_results/CD8-naive/model_results_cells_10_hpc.csv'
#         ,col.names = c("gene","enhancer1","enhancer2","intercept","beta_estimate","beta_pvalue","bootstrap_pvalue"))
# enhancer2 <- read.csv('/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/results/hpc_model_results/CD8-naive/model_results_cells_01_hpc.csv'
#         ,col.names = c("gene","enhancer1","enhancer2","intercept","beta_estimate","beta_pvalue","bootstrap_pvalue"))
# interaction <- read.csv('/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/results/hpc_model_results/CD8-naive/model_results_combined_test_hpc.csv'
#         ,col.names = c("gene", "enhancers", "intercept", "beta.estimate", "beta.pvalue", "bootstrap.pvalue"))


#process interaction matrix
interaction <- interaction %>% mutate(enhancers = as.character(enhancers))
interaction <- interaction %>% separate(enhancers, into = c("enhancer1", "enhancer2"), sep = "_")
interaction <- interaction %>% select(gene, enhancer1, enhancer2, beta.estimate, intercept, beta.pvalue, bootstrap.pvalue)
interaction[0:5,]
# enhancer1[0:5,]

#merge all three matrices into one
merged_df <- merge(interaction, enhancer1, by = c("enhancer1", "enhancer2",'gene'))
merged_df <- merge(merged_df, enhancer2, by = c("enhancer1", "enhancer2",'gene'))
merged_df[0:5,]

# filter out NA values from model results
epistasis.models <- merged_df
epistasis.models <- epistasis.models[complete.cases(epistasis.models), ]

#set column names for model results
colnames(epistasis.models) <- c('enhancer1','enhancer2','gene','beta.estimate.both','intercept.both','beta.pvalue.both','bootstrap.pvalue.both',
   'intercept.10','beta.estimate.10','beta.pvalue.10','bootstrap.pvalue.10',
   'intercept.01','beta.estimate.01','beta.pvalue.01','bootstrap.pvalue.01'
)
epistasis.models[0:5,]

# compute -log(p-value)
print(nrow(epistasis.models))
epistasis.models <- epistasis.models[complete.cases(epistasis.models), ]
print(nrow(epistasis.models))
epistasis.models$neglog10p <- -1 * log10(as.numeric(epistasis.models$bootstrap.pvalue.both))

# mark whether significant or not
epistasis.models$adj.p <- p.adjust(as.numeric(epistasis.models$bootstrap.pvalue.both), method = 'fdr')
epistasis.models$is.significant <- epistasis.models$adj.p < 0.1
print(sum(epistasis.models$is.significant))
cols <- c("TRUE"= 'red', "FALSE" = 'black')
epistasis.models[0:5,]


#scatterplot of enhancer1 + enhancer 2 vs. enhancer_both
# Calculate the sum of beta.estimate_enhancer1 and beta.estimate_enhancer2
# Determine the range for both axes
epistasis.models$beta_sum_enhancers <- epistasis.models$beta.estimate.10 + epistasis.models$beta.estimate.01
axis_range <- range(c(epistasis.models$beta_sum_enhancers, epistasis.models$beta.estimate.both), na.rm = TRUE)
epistasis.models[0:5,]



#figure 1: scatterplot

scatter_plot <- ggplot(epistasis.models, aes(x = beta_sum_enhancers, y = beta.estimate.both)) +
  geom_point(color = "blue", alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +  # y = x line
  coord_equal() +  # Ensure equal scaling on both axes
  xlim(axis_range) +  # Set the x-axis limits
  ylim(axis_range) +  # Set the y-axis limits
  theme_minimal() +
  theme(
    aspect.ratio = 1,  # Make the plot square
    plot.title = element_text(size = 16, face = "bold"),  # Title font size
    axis.title = element_text(size = 14),  # Axis title font size
    axis.text = element_text(size = 12)  # Axis text font size
  ) +
  labs(
    title = "Scatterplot of Enhancer Sum vs Both (Equal Scales, Square, y = x)",
    x = "Sum of Beta Estimates (Enhancer 1 + Enhancer 2)",
    y = "Beta Estimate (Both Enhancers)"
  )

print(scatter_plot)

# Save the plot as a square
ggsave("results/plot/scatterplot_enhancer_sum_vs_both_square.png", scatter_plot, width = 6, height = 6)





#figure 2:box plot 
# Combine the data into a data frame
boxplot_data <- data.frame(
  Model = rep(c("beta_sum_enhancers", "beta_estimate_both"), 
              times = c(length(epistasis.models$beta_sum_enhancers), length(epistasis.models$beta.estimate.both))),
  Beta = c(epistasis.models$beta_sum_enhancers, epistasis.models$beta.estimate.both)
)

boxplot <- ggplot(boxplot_data, aes(x = Model, y = Beta, fill = Model)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Comparison of Beta Distributions",
    x = "Models",
    y = "Beta Values"
  ) +
  theme(
    aspect.ratio = 1.5,                     # Make the plot square
    text = element_text(size = 16),       # General text size
    axis.title = element_text(size = 18),# Axis titles (x and y labels)
    axis.text = element_text(size = 16), # Axis tick labels
    legend.title = element_text(size = 16), # Legend title
    legend.text = element_text(size = 16),  # Legend text
    plot.title = element_text(size = 20, face = "bold") # Plot title
  )

# Print the plot
print(boxplot)

# Save the plot
ggsave("results/plot/beta_comparison_boxplot.png", boxplot)









# #-------------------------------------------------------
# # For each of 10,01,11

#volcano plot for each
volcano_data <- epistasis.models %>%
  select(
    gene,
    neglog10p,
    is.significant,
    beta.estimate.10,
    beta.estimate.01,
    beta.estimate.both
  ) %>%
  pivot_longer(
    cols = starts_with("beta.estimate"),
    names_to = "Enhancer",
    values_to = "Beta"
  )


# Create faceted volcano plot with larger font sizes
faceted_volcano_plot <- ggplot(volcano_data, aes(Beta, neglog10p)) +
  geom_point(aes(color = factor(is.significant)), alpha = 0.5) +
  facet_wrap(~ Enhancer, scales = "free_x", labeller = label_parsed) +
  scale_color_manual(values = cols) +
  theme_classic() +
  labs(
    title = "Volcano Plot Comparison",
    x = "Beta Estimate",
    y = "-log10(P-value)",
    color = "Significant"
  ) +
  theme(
    aspect.ratio = 4,  # Make the plot square
    text = element_text(size = 16),        # General text size
    axis.title = element_text(size = 18), # Axis titles
    axis.text = element_text(size = 14),  # Axis labels
    strip.text = element_text(size = 16), # Facet labels
    plot.title = element_text(size = 20, face = "bold") # Plot title
  )

print(faceted_volcano_plot)
ggsave("results/plot/each/faceted_volcano_plot.png", faceted_volcano_plot)




#Boxplot
boxplot <- ggplot(volcano_data, aes(Enhancer, Beta)) +
  geom_boxplot(aes(fill = factor(is.significant)), alpha = 0.7) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  labs(
    title = "Beta Estimate Distributions",
    x = "Enhancer Type",
    y = "Beta Estimate",
    fill = "Significant"
  ) +
  theme(
    aspect.ratio = 1,  # Make the plot square
    text = element_text(size = 16),        # General text size
    axis.title = element_text(size = 18), # Axis titles
    axis.text = element_text(size = 16),  # Axis labels
    strip.text = element_text(size = 16), # Facet labels
    plot.title = element_text(size = 20, face = "bold") # Plot title
  )
print(boxplot)
ggsave("results/plot/each/boxplot_beta_estimates.png", boxplot)







# #-------------------------------------------------------

# read in RNA-seq matrix, ATAC-seq matrix, and cell-level metadata
rna <- readRDS('/Users/lyricli/Documents/Capstone/Data/PBMC/rna_matrix.rds')
atac <- readRDS('/Users/lyricli/Documents/Capstone/Data/PBMC/atac_matrix.rds')
metadata <- readRDS('/Users/lyricli/Documents/Capstone/Data/PBMC/metafile_raw.rds')
# desired_celltypes <-  config['celltype']

# read in enhancer pairs determined from SCENT
enhancer.pairs <- read.csv('/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/results/enhancer_pairs/CD14-mono/enhancer_pairs_1.csv')

# rna <- rna[, metadata$celltype %in% desired_celltypes]
# atac <- atac[, metadata$celltype %in% desired_celltypes]
# metadata <- metadata[metadata$celltype %in% desired_celltypes,]

# Test with one pair (e.g., i = 1)
i <- 1
enhancer.1 <- enhancer.pairs$enhancer_1[i]
enhancer.2 <- enhancer.pairs$enhancer_2[i]
gene <- enhancer.pairs$gene[i]
gene


# Create categories and data
enhancer_1_activity <- atac[enhancer.1, ]
enhancer_2_activity <- atac[enhancer.2, ]
categories <- ifelse(enhancer_1_activity == 0 & enhancer_2_activity == 0, "00",
              ifelse(enhancer_1_activity == 1 & enhancer_2_activity == 0, "10",
              ifelse(enhancer_1_activity == 0 & enhancer_2_activity == 1, "01", "11")))
plot_data <- data.frame(
  gene_expression = rna[gene, ],
  enhancer_1 = enhancer.1,
  enhancer_2 = enhancer.2,
  category = categories
)
plot_data$normalized_expression <- log1p(plot_data$gene_expression)

# Plot
violin_plot <- ggplot(plot_data, aes(x = category, y = normalized_expression, fill = category)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = paste("Gene Expression of", gene, "\nfor Pair", enhancer.1, " & ",enhancer.2 ),
       x = "Category", y = "Normalized Expression") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(
    aspect.ratio = 1,  # Make the plot square
    text = element_text(size = 16),        # General text size
    axis.title = element_text(size = 18), # Axis titles
    axis.text = element_text(size = 14),  # Axis labels
    strip.text = element_text(size = 16), # Facet labels
    plot.title = element_text(size = 20, face = "bold") # Plot title
  )

print(violin_plot)
ggsave("results/plot/violinplot_enhancer_pair.png", violin_plot)







# #-------------------------------------------------------
# # QQplot for SCENT output 
library(qqman)

# scent <- read.table("/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/results/SCENT_outputs/SCENT_output_1.csv",
#                     header = TRUE, sep = " ")

folder_path <- "/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/results/SCENT_outputs/CD4-mono/"
file_list <- list.files(path = folder_path, pattern = "SCENT_output_.*\\.csv", full.names = TRUE)
read_custom_table <- function(file) {
  read.table(file, header = TRUE, sep = " ")
}
scent <- do.call(rbind, lapply(file_list, read_custom_table))
head(scent)

scent[0:5,]
p_values <- scent$"boot_basic_p"  # or data$boot_basic_p
class(p_values)
# QQ plot using qqman

# Function to generate the QQ plot
generate_qq_plot <- function() {
  par(cex.axis = 1.5, cex.lab = 2)  # Set axis and label sizes
  qq(p_values, main = "", xlab = "", ylab = "")  # Create plot
  abline(0, 1, col = "red", lwd = 2.5)  # Add a thicker red line
  title(main = "QQ Plot for SCENT output: bootstrap p_value", cex.main = 2.5)
  mtext(side = 1, line = 2.5, "Expected -log10(p)", cex = 2)  # X-axis label
  mtext(side = 2, line = 2.5, "Observed -log10(p)", cex = 2)  # Y-axis label
}

# Display the plot on the screen
generate_qq_plot()

# Save the plot to a file
png("results/plot/QQ_plot_SCENT_output.png", width = 700, height = 800, res = 50)
generate_qq_plot()
dev.off()