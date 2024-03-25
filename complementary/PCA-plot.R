#### PCA PLOT ####

file1 <- "/Users/danischiavinato/Desktop/data-from-server/20240320/protein-features"
file2 <- "/Users/danischiavinato/Desktop/data-from-server/20240320/lncrna-features"
file3 <- "/Users/danischiavinato/Desktop/data-from-server/20240320/short-features"

matrix_protein <- data.frame(Dataset = "protein-coding", read.csv(file1, header = TRUE))
matrix_lncrna <- data.frame(Dataset = "lncRNA", read.csv(file2, header = TRUE))
matrix_short_ncrna <- data.frame(Dataset = "short-ncRNA", read.csv(file3, header = TRUE))

all_matrix <- rbind(matrix_protein, matrix_lncrna, matrix_short_ncrna)
all_matrix <- na.omit(all_matrix) 

# Remove the 'Functional' column and store it separately
functional_column <- all_matrix$Functional
dataset_column <- all_matrix$Dataset
feature_matrix <- all_matrix[, -which(names(all_matrix) == "Functional")]
feature_matrix_numeric <- feature_matrix[, sapply(feature_matrix, is.numeric)]

# Perform PCA
pca_result <- prcomp(feature_matrix_numeric, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pc_scores <- as.data.frame(pca_result$x[, 1:2])
pca_df <- cbind(pc_scores, all_matrix[, !names(all_matrix) %in% c("Functional", "Dataset")])

# Add 'Functional' and 'Dataset' columns back to the PCA results data frame
pca_df$Functional <- functional_column
pca_df$Dataset <- dataset_column
pca_df$ColorVariable <- paste(pca_df$Functional, pca_df$Dataset, sep = "-")
functional_colors <- c("Yes" = "blue", "No" = "red")


color_palette <- c("Yes-protein-coding" = "red", "Yes-lncRNA" = "darkorange", "Yes-short-ncRNA" = "pink2",
                   "No-protein-coding" = "blue", "No-lncRNA" = "purple", "No-short-ncRNA" = "deepskyblue")

# Plot 
ggplot(pca_df, aes(x = PC1, y = PC2, color = ColorVariable)) +
    geom_point(size = 1.5, fill = 'white') +
    scale_color_manual(values = color_palette) +  # Specify colors for Functional
    labs(title = "Principal Component Analysis", x = "PC1", y = "PC2", color = "Functional", shape = "Dataset") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 10),    # Increase axis label size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 20),   # Increase legend text size
    legend.key.size = unit(4, "lines")       # Increase legend key size
  )
  



