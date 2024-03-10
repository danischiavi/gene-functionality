# file1 <- "/Users/danischiavinato/Desktop/features/protein-exon2-spearman-matrix"
# file2 <- "/Users/danischiavinato/Desktop/features/protein-exon3-spearman-matrix"
# file3 <- "/Users/danischiavinato/Desktop/features/lncrna-exon1-spearman-matrix"
# file4 <- "/Users/danischiavinato/Desktop/features/lncrna-exon2-spearman-matrix"
# file5 <- "/Users/danischiavinato/Desktop/features/short-ncrna-spearman-matrix"

install.packages("ggplot2")
library(ggplot2)

cor_matrix_protein_exon2 <- read.csv(file1, header = TRUE)
cor_matrix_protein_exon3 <- read.csv(file2, header = TRUE)
cor_matrix_lncrna_exon1 <- read.csv(file3, header = TRUE)
cor_matrix_lncrna_exon2 <- read.csv(file4, header = TRUE)
cor_matrix_short_ncrna <- read.csv(file5, header = TRUE)


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
cor_data_protein_exon2 <- extract_cor_data(cor_matrix_protein_exon2, "Protein_Exon2")
cor_data_protein_exon3 <- extract_cor_data(cor_matrix_protein_exon3, "Protein_Exon3")
cor_data_lncrna_exon1 <- extract_cor_data(cor_matrix_lncrna_exon1, "LncRNA_Exon1")
cor_data_lncrna_exon2 <- extract_cor_data(cor_matrix_lncrna_exon2, "LncRNA_Exon2")
cor_data_short_ncrna <- extract_cor_data(cor_matrix_short_ncrna, "Short_ncRNA")

# Combine all correlation data into a single data frame
all_cor_data <- rbind(cor_data_protein_exon2, cor_data_protein_exon3, 
                      cor_data_lncrna_exon1, cor_data_lncrna_exon2, 
                      cor_data_short_ncrna)


# Create the plot
cor_plot <- ggplot(all_cor_data, aes(x = Feature, y = Correlation, color = Dataset)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.9), size = 0.6) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), position = position_dodge(width = 0.9), width = 0.3) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.25, linetype = "dashed", color = "black") +
  labs(title = "Spearman Correlation Values with Confidence Intervals",
       x = "Feature",
       y = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(cor_plot)