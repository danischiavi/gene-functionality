# file <- "/Users/danischiavinato/Desktop/data-from-server/20240320/short-features"

#### SPEARMAN CORRELATION MATRIX: FUNCTION-FEATURE for each dataset individually #### 
options(repos = "https://cloud.r-project.org")
install.packages("RVAideMemoire")
install.packages("randomForest")

library(tools)
library(RVAideMemoire)
library(randomForest)

# Function declaration to calculate Rho spearman and its confidence intervals 
spearman_correlation <- function(file){

  data <- read.csv(file, header=TRUE)
  
  # Convert Functional "Yes" "No" into numeric values 
  data$Functional <- factor(data$Functional, levels=c("No","Yes")) 
  functional_numeric <- as.numeric(factor(data$Functional)) - 1 

  # Approx NA values (CHECK method) 
  data <- na.roughfix(data)

  # 2. Extract the other numeric columns
  numeric_data <- data[, !names(data) %in% c("Functional")]

  # Calculate the Spearman correlation between "functional" column and each numeric column
  correlations <- sapply(numeric_data, function(col) {
    spearman <- spearman.ci(functional_numeric, col, conf.level = 0.95, nrep = 1000)
    rho <- spearman$estimate[[1]]
    ci_inf <- spearman$conf.int[[1]] 
    ci_sup <- spearman$conf.int[[2]]
    
    correlation_test <- cor.test(functional_numeric, col, method = "spearman")
    p_value <- correlation_test$p.value
    
    return(c(cor = rho, ci_inf = ci_inf, ci_sup = ci_sup, p_value = p_value))
  })
  
  # Return correlation coefficients and confidence intervals as a named list with variable name as the name
  return(correlations)
}

# Initialize empty list to store correlation vectors
correlation_list <- list()

correlation_list[[file]] <- spearman_correlation(file)
  
correlation_matrix_combined <- do.call(rbind, correlation_list)

# Format output file 
filename <- sub("features\\$", "", file)
output_file <- paste0(filename,"-spearman-matrix")
write.csv(correlation_matrix_combined, file = output_file, row.names = TRUE)


### PLOT (all dataset in the same graph) ###

file1 <- "/Users/danischiavinato/Desktop/data-from-server/20240320/protein-features-spearman-matrix"
file2 <- "/Users/danischiavinato/Desktop/data-from-server/20240320/lncrna-features-spearman-matrix"
file3 <- "/Users/danischiavinato/Desktop/data-from-server/20240320/short-features-spearman-matrix"

install.packages("ggplot2")
library(ggplot2)

cor_matrix_protein <- read.csv(file1, header = TRUE)
cor_matrix_lncrna <- read.csv(file2, header = TRUE)
cor_matrix_short_ncrna <- read.csv(file3, header = TRUE)

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


