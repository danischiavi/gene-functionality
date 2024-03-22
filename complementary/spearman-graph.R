file1 <- "/Users/danischiavinato/Desktop/data-from-server/20240320/protein-features-spearman-matrix"
file2 <- "/Users/danischiavinato/Desktop/data-from-server/20240320/lncrna-features-spearman-matrix"
file3 <- "/Users/danischiavinato/Desktop/data-from-server/20240320/short-features-spearman-matrix"

file1 <- "/Users/danischiavinato/Desktop/data-from-server/20240312/protein-exon2-conservation-test.csvspearman-matrix"
file2 <- "/Users/danischiavinato/Desktop/data-from-server/20240312/protein-exon3-conservation-test.csvspearman-matrix"
file3 <- "/Users/danischiavinato/Desktop/data-from-server/20240312/lncrna-exon1-conservation-test.csvspearman-matrix"
file4 <- "/Users/danischiavinato/Desktop/data-from-server/20240312/lncrna-exon2-conservation-test.csvspearman-matrix"
file5 <- "/Users/danischiavinato/Desktop/data-from-server/20240312/short-ncrna-conservation-test.csvspearman-matrix"

install.packages("ggplot2")
library(ggplot2)

cor_matrix_protein <- read.csv(file1, header = TRUE)
cor_matrix_lncrna <- read.csv(file2, header = TRUE)
cor_matrix_short_ncrna <- read.csv(file3, header = TRUE)


# GRAPH 
# Extract correlation values and confidence intervals
extract_cor_data <- function(cor_matrix, dataset_name) {
  cor_values <- cor_matrix[1, -1]
  ci_inf <- cor_matrix[2, -1]
  ci_sup <- cor_matrix[3, -1]
  
  cor_data <- data.frame(
    Dataset = dataset_name,
    Feature = names(cor_values),
    Correlation = as.numeric(as.vector(cor_values)),
    Lower_CI = as.numeric(as.vector(ci_inf)),
    Upper_CI = as.numeric(as.vector(ci_sup))
  )
  
  return(cor_data)
}

# Extract correlation data for each dataset
cor_data_protein <- extract_cor_data(cor_matrix_protein, "Protein-coding")
cor_data_lncrna <- extract_cor_data(cor_matrix_lncrna, "lncRNA")
cor_data_short_ncrna <- extract_cor_data(cor_matrix_short_ncrna, "short-ncRNA")

# Combine all correlation data into a single data frame
all_cor_data <- rbind(cor_data_protein, cor_data_lncrna, cor_data_short_ncrna)

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

# Apply the function to assign groups to each feature
all_cor_data$Group <- sapply(all_cor_data$Feature, assign_group)

# Create the plot
cor_plot <- ggplot(all_cor_data, aes(x = factor(Feature, levels = unique(all_cor_data$Feature)), y = Correlation, color = Group, shape = Dataset)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), position = position_dodge(width = 0.9), width = 0.3) +
  #geom_hline(yintercept = 0.25, linetype = "dashed", color = "black") +
  #geom_hline(yintercept = -0.25, linetype = "dashed", color = "black") +
  labs(title = "Spearman Correlation Features-Function",
       x = "Feature",
       y = "Correlation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 10),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 14),    # Increase axis label size
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 16),   # Increase legend text size
    legend.key.size = unit(2, "lines")       # Increase legend key size
    )

print(cor_plot)

## PCA PLOT ##

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

# Plot PCA
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
  




## PhyloP ##
file <- "/Users/danischiavinato/Desktop/data-from-server/conservation/spearman-corr"

data <- read.csv(file, header = TRUE)

ggplot(data, aes(x = score, y = rho, color = rna)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_inf, ymax = ci_sup), width = 0.1) +
  labs(title = "Correlation Analysis",
       x = "Score Type",
       y = "Correlation Coefficient") +
  theme_minimal()



### FEATURES CORRELATION MATRIX ###
file <- "/Users/danischiavinato/Desktop/data-from-server/20240320/short-features"

data <- read.csv(file, header = TRUE)

data$Functional <- factor(data$Functional, levels=c("No","Yes")) 
functional_numeric <- as.numeric(factor(data$Functional)) - 1 
data$Functional <- functional_numeric
data <- na.roughfix(data)

cor_matrix <- cor(data, method = "spearman")

# Heatmap
library(ggplot2)
library(reshape2)

# Melt the correlation matrix into long format
cor_melted <- melt(cor_matrix)

ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Spearman Correlation") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    plot.title = element_text(size = 12),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 20),    # Increase axis label size
    legend.title = element_text(size = 10),  # Increase legend title size
    legend.text = element_text(size = 10),   # Increase legend text size
    legend.key.size = unit(2, "lines")       # Increase legend key size
  ) +
  labs(title = "Features Spearman Correlation - short-ncRNA", x = "Features", y = "Features")