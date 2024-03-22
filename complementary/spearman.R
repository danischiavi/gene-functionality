# file <- "/Users/danischiavinato/Desktop/data-from-server/20240320/short-features"



options(repos = "https://cloud.r-project.org")
install.packages("RVAideMemoire")
install.packages("randomForest")
#install.packages("tools", lib="packages")
library(tools)
library(RVAideMemoire)
library(randomForest)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

# Read input file from Snakemake input file directive
# file <- input$file

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



### RANDOM NUMBER - NEUTRAL PREDICTOR ###
i=0
while [ $i -lt 885 ]; do 
  echo $RANDOM >> random_numbers.tmp
  (( i++))
done 

# Paste the generated random numbers as a new column into the existing file
paste -d ',' random_numbers.tmp $file > lncrna-features

# Clean up temporary file
rm random_numbers.tmp
