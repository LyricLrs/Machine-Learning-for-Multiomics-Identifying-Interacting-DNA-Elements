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

# read in epistasis models from batch
enhancer1 <- read.csv('/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/workflow/scripts/model_results/model_results_cells_11 10.27.44 AM.csv'
        ,col.names = c("gene", "enhancer1", "intercept", "beta.estimate", "beta.pvalue", "bootstrap.pvalue"))
enhancer2 <- read.csv('/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/workflow/scripts/model_results/model_results_cells_21 10.27.45 AM.csv'
        ,col.names = c("gene", "enhancer2", "intercept", "beta.estimate", "beta.pvalue", "bootstrap.pvalue"))
interaction <- read.csv('/Users/lyricli/Documents/Capstone/CDS-2024-Fall-Capstone/workflow/scripts/model_results/model_results_combined.csv'
        ,col.names = c("gene", "enhancers", "intercept", "beta.estimate", "beta.pvalue", "bootstrap.pvalue"))


#process interaction matrix
interaction <- interaction %>% mutate(enhancers = as.character(enhancers))
interaction <- interaction %>% separate(enhancers, into = c("enhancer1", "enhancer2"), sep = "_")
interaction <- interaction %>% select(gene, enhancer1, enhancer2, beta.estimate, intercept, beta.pvalue, bootstrap.pvalue)

#merge all three matrices into one
merged_df <- enhancer1 %>% inner_join(interaction, by = "enhancer1", suffix = c("_enhancer1", "_both"))
colnames(enhancer2) <- ifelse(
  colnames(enhancer2) %in% c("gene", "enhancer2"),
  colnames(enhancer2),
  paste0(colnames(enhancer2), "_enhancer2")
)
merged_df <- merged_df %>% inner_join(enhancer2, by = "enhancer2")
merged_df[0:5,]


# filter out NA values from model results
epistasis.models <- merged_df
epistasis.models <- epistasis.models[complete.cases(epistasis.models), ]
epistasis.models[0:5,]
# set column names for model results
# colnames(epistasis.models) <- c(
#     "gene", "enhancer1", "intercept", "beta.estimate", "beta.pvalue", "bootstrap.pvalue"
# )

# compute -log(p-value)
print(nrow(epistasis.models))
epistasis.models <- epistasis.models[complete.cases(epistasis.models), ]
print(nrow(epistasis.models))
epistasis.models$neglog10p <- -1 * log10(as.numeric(epistasis.models$bootstrap.pvalue_both))

# mark whether significant or not
epistasis.models$adj.p <- p.adjust(as.numeric(epistasis.models$bootstrap.pvalue_both), method = 'fdr')
epistasis.models$is.significant <- epistasis.models$adj.p < 0.1
print(sum(epistasis.models$is.significant))
cols <- c("TRUE"= 'red', "FALSE" = 'black')



#scatterplot of enhancer1 + enhancer 2 vs. enhancer_both
# Calculate the sum of beta.estimate_enhancer1 and beta.estimate_enhancer2
# Determine the range for both axes
axis_range <- range(c(epistasis.models$beta_sum_enhancers, epistasis.models$beta.estimate_both), na.rm = TRUE)

# Create the scatterplot with fixed equal scales and axis limits
scatter_plot <- ggplot(epistasis.models, aes(x = beta_sum_enhancers, y = beta.estimate_both)) +
  geom_point(color = "blue", alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +  # y = x line
  coord_equal() +  # Ensure equal scaling on both axes
  xlim(axis_range) +  # Set the x-axis limits
  ylim(axis_range) +  # Set the y-axis limits
  theme_minimal() +
  theme(aspect.ratio = 1) +  # Make the plot square
  labs(
    title = "Scatterplot of Enhancer Sum vs Both (Equal Scales, Square, y = x)",
    x = "Sum of Beta Estimates (Enhancer 1 + Enhancer 2)",
    y = "Beta Estimate (Both Enhancers)"
  )

# Save the plot as a square
ggsave("CDS-2024-Fall-Capstone/results/plot/scatterplot_enhancer_sum_vs_both_square.png", scatter_plot, width = 6, height = 6)

# Display the plot
print(scatter_plot)

# #volcano plot for each
# volcano_data <- epistasis.models %>%
#   select(
#     gene_both,
#     neglog10p,
#     is.significant,
#     beta.estimate_enhancer1,
#     beta.estimate_enhancer2,
#     beta.estimate_both
#   ) %>%
#   pivot_longer(
#     cols = starts_with("beta.estimate"),
#     names_to = "Enhancer",
#     values_to = "Beta"
#   )

# # Create faceted volcano plot
# ggplot(volcano_data, aes(Beta, neglog10p)) +
#   geom_point(aes(color = factor(is.significant)), alpha = 0.5) +
#   facet_wrap(~ Enhancer, scales = "free_x", labeller = label_parsed) +
#   scale_color_manual(values = cols) +
#   theme_classic() +
#   labs(
#     title = "Volcano Plot Comparison",
#     x = "Beta Estimate",
#     y = "-log10(P-value)",
#     color = "Significant"
#   )

# ggsave("results/plot/faceted_volcano_plot.png", plot)



# #Beta Estimates Across Genes
# # Melt the merged_df to long format
# line_plot_df <- melt(
#   epistasis.models,
#   id.vars = "gene_both",  # Keep 'gene_both' as the identifier
#   measure.vars = c("beta.estimate_enhancer1", "beta.estimate_enhancer2", "beta.estimate_both"),  # Columns to melt
#   variable.name = "Type",  # Name for the new 'Type' column
#   value.name = "Beta Estimate"  # Name for the new 'Beta Estimate' column
# )

# # Plot the data
# plot <- ggplot(line_plot_df, aes(x = gene_both, y = `Beta Estimate`, color = Type, group = Type)) +
#   geom_line() +                # Add lines
#   geom_point() +               # Add points
#   theme_minimal() +            # Clean theme
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
#   labs(
#     title = "Beta Estimates Across Genes",
#     x = "Gene",
#     y = "Beta Estimate"
#   )

# ggsave("results/plot/Beta_estimates_across_genes.png", plot)



# #Boxplot
# plot <- ggplot(volcano_data, aes(Enhancer, Beta)) +
#   geom_boxplot(aes(fill = factor(is.significant)), alpha = 0.7) +
#   scale_fill_manual(values = cols) +
#   theme_classic() +
#   labs(
#     title = "Beta Estimate Distributions",
#     x = "Enhancer Type",
#     y = "Beta Estimate",
#     fill = "Significant"
#   )

# ggsave("results/plot/boxplot_beta_estimates.png", plot)





# #-------------------------------------------------------
# # Compute pairwise differences
# contrast_data <- epistasis.models %>%
#   mutate(
#     diff_1_both = beta.estimate_enhancer1 - beta.estimate_both,
#     diff_2_both = beta.estimate_enhancer2 - beta.estimate_both
#   )

# # Prepare data for comparison plot
# contrast_data <- contrast_data %>%
#   select(gene_both, diff_1_both, diff_2_both) %>%
#   pivot_longer(
#     cols = starts_with("diff"),
#     names_to = "Comparison",
#     values_to = "Difference"
#   )

# # Plot pairwise differences
# plot <- ggplot(contrast_data, aes(x = Comparison, y = Difference, fill = Difference)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_minimal() +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   labs(
#     title = "Pairwise Beta Estimate Differences",
#     x = "Comparison",
#     y = "Difference",
#     fill = "Magnitude"
#   )

# ggsave("results/plot/pairwise_differences.png", plot)


# # Prepare data for box plot
# boxplot_data <- epistasis.models %>%
#   mutate(
#     diff_1_both = beta.estimate_enhancer1 - beta.estimate_both,
#     diff_2_both = beta.estimate_enhancer2 - beta.estimate_both
#   ) %>%
#   select(gene_both, diff_1_both, diff_2_both) %>%
#   pivot_longer(
#     cols = starts_with("diff"),
#     names_to = "Comparison",
#     values_to = "Difference"
#   )

# # Create box plot
# boxplot <- ggplot(boxplot_data, aes(x = Comparison, y = Difference)) +
#   geom_boxplot(fill = "blue", alpha = 0.7, outlier.color = "red") +
#   theme_minimal() +
#   labs(
#     title = "Box Plot of Pairwise Beta Estimate Differences",
#     x = "Comparison",
#     y = "Difference"
#   )

# # Save the plot
# ggsave("results/plot/boxplot_pairwise_differences.png", boxplot)


