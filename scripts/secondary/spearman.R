# file <- "/Users/danischiavinato/Desktop/features/protein-exon2-features.csv"
# file <- "/Users/danischiavinato/Desktop/features/protein-exon3-features.csv"
# file <- "/Users/danischiavinato/Desktop/features/lncrna-exon1-features.csv"

options(repos = "https://cloud.r-project.org")
install.packages("RVAideMemoire", lib="packages")
library(RVAideMemoire)
install.packages("randomForest", lib="packages")
library(randomForest)
#install.packages("tools", lib="packages")
library(tools)


args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

# Read input file from Snakemake input file directive
# file <- input$file

# Format output file 
filename <- sub("features\\.csv$", "", file)
output_file <- paste0(filename,"spearman-matrix")


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
    return(c( cor = rho, ci_inf = ci_inf, ci_sup = ci_sup ))
  })
  
  # Return correlation coefficients and confidence intervals as a named list with variable name as the name
  return(correlations)
}


# Initialize empty list to store correlation vectors
correlation_list <- list()

correlation_list[[file]] <- spearman_correlation(file)
  
correlation_matrix_combined <- do.call(rbind, correlation_list)

# Write the correlation matrix to a CSV file
write.csv(correlation_matrix_combined, file = output_file, row.names = TRUE)

