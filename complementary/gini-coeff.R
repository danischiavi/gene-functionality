gini_matrix <- read.csv("/Users/danischiavinato/Desktop/data-from-server/20240320/protein-features-importance", header = TRUE)
gini_matrix <- read.csv("/Users/danischiavinato/Desktop/data-from-server/20240320/lncrna-features-importance", header = TRUE)
gini_matrix <- read.csv("/Users/danischiavinato/Desktop/data-from-server/20240320/short-features-importance", header = TRUE)

mean_gini <- rowMeans(gini_matrix[, -1])
num_models <- ncol(gini_matrix) - 1  # Number of models
sd_gini <- apply(gini_matrix[, -1], 1, sd)  # Standard deviation of Gini coefficients for each feature
sem_gini <- sd_gini / sqrt(num_models)
# Create a data frame for plotting

importance_df <- data.frame(Feature = gini_matrix[, 1], MeanGini = mean_gini, SEM = sem_gini)

# Sort the data frame by importance
importance_df <- importance_df[order(importance_df$MeanGini, decreasing = TRUE), ]

# custom_order <- c("X241w_PP_mean", "X100w_PP_mean", "X241w_PP_max", "X100w_PP_max") --> to use as level 
random <- c("Neutral")
intrinsic <- c("GC.")
conservation <- c("phyloP_mean", "phyloP_max", "phastCons_mean", "phastCons_max", "GERP_mean", "GERP_max")
transcriptome <- c("RPKM_tissue", "MRD_tissue", "RPKM_primary.cell", "MRD_primary.cell")
specific <- c("Interaction_min", "Interaction_ave", "Fickett_score", "RNAcode_score", "RNAalifold_score", "Accessibility", "Max_covariance", "Min_covariance_Eval", "MFE")
repeats <- c("copy_number", "Dfam_min", "Dfam_sum")
population <- c("SNP_density", "MAF_avg")

# Combine groups into a list
groups <- list(intrinsic = intrinsic, conservation = conservation, transcriptome = transcriptome, specific = specific, repeats = repeats, population = population, neutral = random)
assign_group <- function(feature_name) {
  for (group_name in names(groups)) {
    if (feature_name %in% groups[[group_name]]) {
      return(group_name)
    }
  }
  return("Other")
}
importance_df$Group <- sapply(importance_df$Feature, assign_group)

# Plot Gini coefficient importance
library(ggplot2)
ggplot(importance_df, aes(x = reorder(Feature, MeanGini), y = MeanGini, color = Group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = MeanGini - SEM, ymax = MeanGini + SEM), width = 0.2, color = "black") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "pink") + #position of the neutral in y axis from the bottom 
  coord_flip() +
  labs(title = "Gini coefficient importance - short-ncRNA",
       x = "Feature",
       y = "Mean Gini Coefficient over 100 runs Importance") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 15),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 15),    # Increase axis label size
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 15),   # Increase legend text size
    legend.key.size = unit(4, "lines")       # Increase legend key size
  )