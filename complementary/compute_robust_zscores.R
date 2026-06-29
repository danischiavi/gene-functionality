###############
# Utility script to Compute robust z-scores for all numeric features separated by Gene Type.
###############
#options(repos = c(CRAN = "https://cran.r-project.org"))
library(dplyr)
library(ggplot2)
library(tidyr)
#install.packages("patchwork", dependencies = TRUE)
library(patchwork)
library(scales)

source("complementary/load_gene_functionality_features.R")
feature_matrix <- load_gene_functionality_features()

# Load epigenetic features (load_epigenetic_features.R not yet present — skipped)
#source("complementary/load_epigenetic_features.R")
#feature_matrix_epigenetic <- load_epigenetic_features()

# Attach RNAID positionally from the gene-functionality feature files (the
# canonical source of the per-row identifier; the histone files behind
# feature_matrix_epigenetic carry no identifier of their own). Both loaders
# rbind the same 10 groups in the same order with no filtering/sampling, so
# row order lines up across feature_matrix, feature_matrix_epigenetic, and
# this file list.
source("complementary/utils.R")
GENE_FUNCTIONALITY_FILES <- c(FUNC_PROT_EXON2_FEATURES_FILE, FUNC_PROT_EXON3_FEATURES_FILE,
                               FUNC_LNCRNA_EXON1_FEATURES_FILE, FUNC_LNCRNA_EXON2_FEATURES_FILE, FUNC_SNCRNA_FEATURES_FILE,
                               NC_PROT_EXON2_FEATURES_FILE, NC_PROT_EXON3_FEATURES_FILE,
                               NC_LNCRNA_EXON1_FEATURES_FILE, NC_LNCRNA_EXON2_FEATURES_FILE, NC_SNCRNA_FEATURES_FILE)
rnaids <- load_rnaids(GENE_FUNCTIONALITY_FILES, "ID")
stopifnot(length(rnaids) == nrow(feature_matrix))
feature_matrix$RNAID <- rnaids
#stopifnot(length(rnaids) == nrow(feature_matrix_epigenetic))
#feature_matrix_epigenetic$RNAID <- rnaids

# Function to compute robust z-scores
compute_robust_zscores_optimized <- function(positive_matrix, negative_matrix, verbose = FALSE) {
  # Use NEGATIVE CONTROLS median absolute deviations (MAD) and medians to compute robust z-scores
  mad_values <- apply(negative_matrix, 2, mad, na.rm = TRUE)  # Calculate MAD for each column in the negative control set
  median_values <- apply(negative_matrix, 2, median, na.rm = TRUE)  # Calculate median for each column in the negative control set
  meanAD_values <- apply(negative_matrix, 2, function(col) mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE))
  
  # Handle case where MAD is 0 or very small by using meanAD instead
  mad_values <- apply(negative_matrix, 2, function(col) {
    mad_value <- mad(col, na.rm = TRUE)
    if (mad_value <= 1e-6) {  # If MAD is 0 or very small, calculate mean absolute deviation (meanAD)
      meanAD <- mean(abs(col - median(col, na.rm = TRUE)), na.rm = TRUE)
      warning(sprintf("Column with very small MAD of %f encountered: Using meanAD value of %f.\n", mad_value, meanAD))
      return(1.2533 * meanAD)  # Use meanAD if it is greater than the threshold
    } else {
      return(mad_value)  # Use MAD if it is sufficient
    }
  })
  
  # Calculate z-scores: subtract the median and divide by MAD for each feature
  zscores <- sweep(positive_matrix, 2, median_values, "-")  # Center each column by subtracting the median
  zscores <- sweep(zscores, 2, mad_values, "/")  # Scale each column by dividing by MAD
  
  print("\n\nMedian values:\n")
  print(median_values)
  print("MAD values:\n")
  print(mad_values)
  print("meanAD values:\n")
  print(meanAD_values)
  
  zscores
}

# Run z-score computation function. Receives a dataframe and an output directory
get_robust_zscores <- function(data, output_dir) {
  # Capture command line arguments
  #args <- commandArgs(trailingOnly = TRUE)
  #if (length(args) < 2) {
  #  stop("Please provide the path to the data file as an argument and output folder.")  # Ensure a file path is provided
  #}
  #names(feature_matrix)[1] <- "row"
  # Load and preprocess data
  #data_file <- args[1]  # Get the file path from the command line arguments
  #data <- read.csv(data_file, header = TRUE, check.names = TRUE)# Load the data from the specified CSV file
  #data <- feature_matrix %>% dplyr::select(-Dataset,-ID,-Functional,-Chromosome,-Start,-End,-Sequence)
  #data <- data_numeric_epigenetic_sample
  
  # Check if subset option is provided
  #subset_value <- as.numeric(gsub("--subset=", "", args[grep("--subset=", args)]))
  #if (is.na(subset_value)) {
  #subset_value <- nrow(data)  # Use nrow to indicate no subsetting by default
  #subset_value <- 100
  #}
  
  # Get output directory
  #output_dir <- args[2]
  #output_dir <- "../results/epigenetic_zscores"
  
  # Convert to numeric
  data_numeric <- data %>%
    dplyr::select(-Dataset, -RNAID) %>%
    sapply(function(feature) as.numeric(as.character(feature))) %>%
    as.data.frame()
  #summary(data_numeric)
  # Restore Dataset and RNAID columns
  data_numeric$Dataset <- data$Dataset
  data_numeric$RNAID <- data$RNAID
  #unique(data_numeric$Dataset)
  #table(data_numeric$Dataset)
  #Protein coding
  # Separate the data into positive and negative feature matrices
  protein_positive_feature_matrix <- data_numeric %>% 
    filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") #%>% 
  #dplyr::select(-Dataset)
  protein_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "protein-exon2-negative-control" | Dataset == "protein-exon3-negative-control")# %>% 
  #dplyr::select(-Dataset)
  
  # Select random samples from positive and negative datasets
  #set.seed(2025)  # Set seed for reproducibility
  #protein_positive_feature_matrix <- protein_positive_feature_matrix %>% sample_n(min(2 * subset_value, nrow(protein_positive_feature_matrix)))  # Subset rows from positive data
  #protein_negative_feature_matrix <- protein_negative_feature_matrix %>% sample_n(min(2 * 10 * subset_value, nrow(protein_negative_feature_matrix)))  # Subset rows from negative data
  
  #print(summary(protein_positive_feature_matrix %>% select(coding_potential, Max_covariance)))
  #print(summary(protein_negative_feature_matrix %>% select(coding_potential, Max_covariance)))  
  #print(apply(protein_positive_feature_matrix, 2, sd, na.rm = TRUE))
  #apply(protein_negative_feature_matrix, 2, sd, na.rm = TRUE)
  #print(sd(protein_negative_feature_matrix$RPKM_tissue, na.rm = TRUE))
  # Compute z-scores
  protein_functional_z_scores <- compute_robust_zscores_optimized(protein_positive_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                                  protein_negative_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                                  verbose)
  protein_negative_z_scores <- compute_robust_zscores_optimized(protein_negative_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                                protein_negative_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                                verbose)
  
  #print(count(protein_positive_feature_matrix))
  #print(summary(protein_negative_feature_matrix))
  #print(sd(protein_functional_z_scores$RPKM_tissue, na.rm = TRUE))
  #print(sd(protein_negative_z_scores$RPKM_tissue, na.rm = TRUE))
  
  #Restore Dataset and RNAID
  protein_functional_z_scores$RNAID <- protein_positive_feature_matrix$RNAID
  protein_functional_z_scores$Dataset <- protein_positive_feature_matrix$Dataset
  protein_functional_z_scores <- protein_functional_z_scores %>% dplyr::relocate(RNAID, Dataset)
  protein_negative_z_scores$RNAID <- protein_negative_feature_matrix$RNAID
  protein_negative_z_scores$Dataset <- protein_negative_feature_matrix$Dataset
  protein_negative_z_scores <- protein_negative_z_scores %>% dplyr::relocate(RNAID, Dataset)

  protein_exon2_functional_z_scores <- protein_functional_z_scores %>% 
                                            filter(Dataset == "protein-coding-exon2")
  protein_exon3_functional_z_scores <- protein_functional_z_scores %>% 
                                                  filter(Dataset == "protein-coding-exon3")

  protein_exon2_negative_z_scores <- protein_negative_z_scores %>% 
                                          filter(Dataset == "protein-exon2-negative-control")
  protein_exon3_negative_z_scores <- protein_negative_z_scores %>% 
                                                filter(Dataset == "protein-exon3-negative-control")
  #Join in one dataframe
  protein_zscores_all <- rbind(protein_functional_z_scores,protein_negative_z_scores)
  # Save z-scores to CSV
  z_scores_file <- file.path(output_dir, "functional-protein-exon2-dataset-zscores.csv")
  write.csv(protein_exon2_functional_z_scores, z_scores_file, row.names = FALSE)
  z_scores_file <- file.path(output_dir, "functional-protein-exon3-dataset-zscores.csv")
  write.csv(protein_exon3_functional_z_scores, z_scores_file, row.names = FALSE)
  
  z_scores_file <- file.path(output_dir, "protein-exon2-negative-control-dataset-zscores.csv")
  write.csv(protein_exon2_negative_z_scores, z_scores_file, row.names = FALSE)
  z_scores_file <- file.path(output_dir, "protein-exon3-negative-control-dataset-zscores.csv")
  write.csv(protein_exon3_negative_z_scores, z_scores_file, row.names = FALSE)
  
  # Save raw data subset
  subset_file <- file.path(output_dir, "mrna_features.csv")
  write.csv(protein_positive_feature_matrix, subset_file, row.names = FALSE)
  
  subset_file <- file.path(output_dir, "mrna_negative_control_features.csv")
  write.csv(protein_negative_feature_matrix, subset_file, row.names = FALSE)
  
  #lncRNA
  # Separate the data into positive and negative feature matrices
  lncrna_positive_feature_matrix <- data_numeric %>% 
    filter(Dataset == "lncrna-exon1" | Dataset == "lncrna-exon2") #%>% 
  #dplyr::select(-Dataset)
  lncrna_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "lncrna-exon1-negative-control" | Dataset == "lncrna-exon2-negative-control") #%>% 
  #dplyr::select(-Dataset)
  
  # Select random samples from positive and negative datasets
  #set.seed(2025)  # Set seed for reproducibility
  #lncrna_positive_feature_matrix <- lncrna_positive_feature_matrix %>% sample_n(min(2 * subset_value, nrow(lncrna_positive_feature_matrix)))  # Subset rows from positive data
  #lncrna_negative_feature_matrix <- lncrna_negative_feature_matrix %>% sample_n(min(2 * 10 * subset_value, nrow(lncrna_negative_feature_matrix)))  # Subset rows from negative data
  
  # Compute z-scores
  lncrna_functional_z_scores <- compute_robust_zscores_optimized(lncrna_positive_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                                 lncrna_negative_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                                 verbose)
  lncrna_negative_z_scores <- compute_robust_zscores_optimized(lncrna_negative_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                               lncrna_negative_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                               verbose)
  
  #print(summary(lncrna_positive_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #count(lncrna_positive_feature_matrix)
  #print(summary(lncrna_negative_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #count(lncrna_negative_feature_matrix)
  
  #Restore Dataset and RNAID
  lncrna_functional_z_scores$RNAID <- lncrna_positive_feature_matrix$RNAID
  lncrna_functional_z_scores$Dataset <- lncrna_positive_feature_matrix$Dataset
  lncrna_functional_z_scores <- lncrna_functional_z_scores %>% dplyr::relocate(RNAID, Dataset)
  lncrna_negative_z_scores$RNAID <- lncrna_negative_feature_matrix$RNAID
  lncrna_negative_z_scores$Dataset <- lncrna_negative_feature_matrix$Dataset
  lncrna_negative_z_scores <- lncrna_negative_z_scores %>% dplyr::relocate(RNAID, Dataset)

  lncrna_exon1_functional_z_scores <- lncrna_functional_z_scores %>% 
                                                 filter(Dataset == "lncrna-exon1")
  lncrna_exon2_functional_z_scores <- lncrna_functional_z_scores %>% 
                                                 filter(Dataset == "lncrna-exon2")
  lncrna_exon1_negative_z_scores <- lncrna_negative_z_scores %>% 
                                               filter(Dataset == "lncrna-exon1-negative-control")
  lncrna_exon2_negative_z_scores <- lncrna_negative_z_scores %>% 
                                               filter(Dataset == "lncrna-exon2-negative-control")
  
  # Save z-scores to CSV
  z_scores_file <- file.path(output_dir, "functional-lncrna-exon1-dataset-zscores.csv")
  write.csv(lncrna_exon1_functional_z_scores, z_scores_file, row.names = FALSE)
  z_scores_file <- file.path(output_dir, "functional-lncrna-exon2-dataset-zscores.csv")
  write.csv(lncrna_exon2_functional_z_scores, z_scores_file, row.names = FALSE)
  
  z_scores_file <- file.path(output_dir, "lncrna-exon1-negative-control-dataset-zscores.csv")
  write.csv(lncrna_exon1_negative_z_scores, z_scores_file, row.names = FALSE)
  z_scores_file <- file.path(output_dir, "lncrna-exon2-negative-control-dataset-zscores.csv")
  write.csv(lncrna_exon2_negative_z_scores, z_scores_file, row.names = FALSE)
  
  # Save raw subset feature to CSV
  subset_file <- file.path(output_dir, "lncrna_features.csv")
  write.csv(lncrna_positive_feature_matrix, subset_file, row.names = FALSE)
  
  subset_file <- file.path(output_dir, "lncrna_negative_control_features.csv")
  write.csv(lncrna_negative_feature_matrix, subset_file, row.names = FALSE)
  
  #sncRNA
  # Separate the data into positive and negative feature matrices
  sncrna_positive_feature_matrix <- data_numeric %>%
    filter(Dataset == "short-ncrna") #%>%
  #dplyr::select(-Dataset)
  sncrna_negative_feature_matrix <- data_numeric %>% 
    filter(Dataset == "short-ncrna-negative-control") #%>% 
  #dplyr::select(-Dataset)
  
  # Select random samples from positive and negative datasets
  #set.seed(2025)  # Set seed for reproducibility
  #sncrna_positive_feature_matrix <- sncrna_positive_feature_matrix %>% sample_n(min(subset_value, nrow(sncrna_positive_feature_matrix)))  # Subset rows from positive data
  #sncrna_negative_feature_matrix <- sncrna_negative_feature_matrix %>% sample_n(min(10 * subset_value, nrow(sncrna_negative_feature_matrix)))  # Subset rows from negative data
  
  # Compute z-scores
  sncrna_functional_z_scores <- compute_robust_zscores_optimized(sncrna_positive_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                                 sncrna_negative_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                                 verbose)
  sncrna_negative_z_scores <- compute_robust_zscores_optimized(sncrna_negative_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                               sncrna_negative_feature_matrix %>% dplyr::select(-Dataset, -RNAID),
                                                               verbose)
  
  #print(summary(sncrna_positive_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #count(sncrna_positive_feature_matrix)
  #print(summary(sncrna_negative_feature_matrix %>% select(H3K9ac_MaxScaledSignal,H3K36me3_MaxScaledSignal,H3K79me2_MaxScaledSignal,chrm_acc_MaxScaledSignal)))
  #count(sncrna_negative_feature_matrix)
  
  #Restore Dataset and RNAID
  sncrna_functional_z_scores$RNAID <- sncrna_positive_feature_matrix$RNAID
  sncrna_functional_z_scores$Dataset <- sncrna_positive_feature_matrix$Dataset
  sncrna_functional_z_scores <- sncrna_functional_z_scores %>% dplyr::relocate(RNAID, Dataset)
  sncrna_negative_z_scores$RNAID <- sncrna_negative_feature_matrix$RNAID
  sncrna_negative_z_scores$Dataset <- sncrna_negative_feature_matrix$Dataset
  sncrna_negative_z_scores <- sncrna_negative_z_scores %>% dplyr::relocate(RNAID, Dataset)

  sncrna_functional_z_scores <- sncrna_functional_z_scores %>% 
                                           filter(Dataset ==  "short-ncrna")
  sncrna_negative_z_scores <- sncrna_negative_z_scores %>% 
                                         filter(Dataset == "short-ncrna-negative-control")
  # Save z-scores to CSV
  z_scores_file <- file.path(output_dir, "functional-short-ncrna-dataset-zscores.csv")
  write.csv(sncrna_functional_z_scores, z_scores_file, row.names = FALSE)
  
  z_scores_file <- file.path(output_dir, "short-ncrna-negative-control-dataset-zscores.csv")
  write.csv(sncrna_negative_z_scores, z_scores_file, row.names = FALSE)
  
  # Save raw subset to CSV
  subset_file <- file.path(output_dir, "sncrna_features.csv")
  write.csv(sncrna_positive_feature_matrix, subset_file, row.names = FALSE)
  
  subset_file <- file.path(output_dir, "sncrna_negative_control_features.csv")
  write.csv(sncrna_negative_feature_matrix, subset_file, row.names = FALSE)
}

# Compute gene functionality z-scores and save to disk:
get_robust_zscores(feature_matrix, ZSCORES_DIR)

# Compute epigenetic z-scores and save to disk:
#get_robust_zscores(feature_matrix_epigenetic, EPIGENETIC_ZSCORES_DIR)
