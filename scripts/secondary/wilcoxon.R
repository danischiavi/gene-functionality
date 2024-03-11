file <- commandArgs(trailingOnly = TRUE)

options(repos = "https://cloud.r-project.org")
install.packages("randomForest")
library(randomForest)
install.packages(tools)
library(tools)

# Format output file 
filename <- sub("features\\.csv$", "", file)
output_file <- paste0(filename,"wilcoxon-matrix")


# Function declaration to calculate Rho spearman and its confidence intervals 
wilcoxon_test <- function(file){
  
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
    wilcoxon <- wilcox.test(functional_numeric, col, conf.level = 0.95, conf.int = TRUE)
    W <- wilcoxon$statistic[[1]]
    p_value <- wilcoxon$p.value[1]
    ci_inf <- wilcoxon$conf.int[[1]] 
    ci_sup <- wilcoxon$conf.int[[2]]
    return(c( cor = W, pvalue = p_value, ci_inf = ci_inf, ci_sup = ci_sup ))
    })
  
  # Return correlation coefficients and confidence intervals as a named list with variable name as the name
  return(correlations)
}


# Initialize empty list to store correlation vectors
correlation_list <- list()

correlation_list[[file]] <- wilcoxon_test(file)

correlation_matrix_combined <- do.call(rbind, correlation_list)

# Write the correlation matrix to a CSV file
write.csv(correlation_matrix_combined, file = output_file, row.names = TRUE)


