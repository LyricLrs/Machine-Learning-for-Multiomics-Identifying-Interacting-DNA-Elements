directory_path <- "/Users/huajingru/Desktop/Fall_2024/Capstone/CD14-Mono/epistasis_models"

# List all files that match the pattern
model_files <- list.files(directory_path, pattern = "11_00_\\d+\\.csv$", full.names = TRUE)

# Read all of them into a list
model_data_list <- lapply(model_files, read.csv)

# Combine them into a single data frame (assuming they have the same columns)
model_results <- do.call(rbind, model_data_list)

# # Separate the second column into 'enhancer_1' and 'enhancer_2'
# model_results <- model_results %>%
#   tidyr::separate(
#     col = 2, # The second column to split
#     into = c("enhancer_1", "enhancer_2"), # New column names
#     sep = "_" # Separator between enhancer_1 and enhancer_2
#   )
# 
# # Rename all other columns
# colnames(model_results) <- c(
#   "gene",          # First column
#   "enhancer_1",    # Newly created column
#   "enhancer_2",    # Newly created column
#   "intercept",     # Third column renamed
#   "beta_estimate", # Fourth column renamed
#   "beta_pvalue",   # Fifth column renamed
#   "bootstrap_pvalue" # Sixth column renamed
# )

zero_expression_summary <- read.csv("per_pair_population_summary.csv")

merged_data <- merge(
  model_results,
  zero_expression_summary,
  by.x = c("enhancer_1", "enhancer_2", "gene"),
  by.y = c("enhancer_1", "enhancer_2", "gene")
)


large_estimates <- merged_data %>%
  filter(abs(beta_estimate) > 1)

if (nrow(large_estimates) > 0) {
  large_estimates_analysis <- large_estimates %>%
    mutate(
      all_zero_expression = total_cells > 0 & zero_expression_count == total_cells,
      zero_expression_fraction = zero_expression_percentages 
    )
  
  # Print the relevant columns
  print(large_estimates_analysis[, c(
    "enhancer_1", "enhancer_2", "gene",
    "beta_estimate", "total_cells",
    "zero_expression_percentages", "all_zero_expression"
  )])
  }
  
  # Save the analysis to a file
write.csv(large_estimates_analysis, "/Users/huajingru/Desktop/Fall_2024/Capstone/CD14-Mono/epistasis_models_11_00_analysis.csv", row.names = FALSE)
