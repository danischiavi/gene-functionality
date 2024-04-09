# Load necessary library
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
library(boot) # Used to obtain confidence intervals


###################################################################
# Sample data vectors for positive and negative controls LONG NCRNA
lncrna_file <- "../data/histone_feature/H3k36me3/H3k36me3_lncrna_matrix.csv"
data_lncrna <- read.csv(lncrna_file, header=TRUE, sep = "\t")

# Calculate basic statistics
summary(data_lncrna$AvgSignal)
# Statistics split by 'Functional' category
by(lncrna_file1$AvgSignal, data_lncrna$Functional, summary) 

# Split data by 'Functional' for t-test
negative_control <- data_lncrna[data_lncrna$Functional == "No", "AvgSignal"]
positive_data <- data_lncrna[data_lncrna$Functional == "Yes", "AvgSignal"]

# Kolmogorov-Smirnov Test (Normality)
ks.test(negative_control, "pnorm", mean(negative_control), sd(negative_control))
ks.test(positive_data, "pnorm", mean(positive_data), sd(positive_data))

# Student's t-test (assuming equal variances)
#t.test(negative_control, positive_data, var.equal = TRUE)

# Mann-Whitney U Test (aka Wilcoxon Rank-Sum Test)
wilcox.test(negative_control, positive_data)
############################################



####################################################################
# Sample data vectors for positive and negative controls SHORT NCRNA
sncrna_file <- "../data/histone_feature/H3k36me3/H3k36me3_short_ncrna_matrix.csv"
data_sncrna <- read.csv(sncrna_file, header=TRUE, sep="\t")

# Calculate basic statistics
summary(data_sncrna$AvgSignal)
# Statistics split by 'Functional' category
by(data_sncrna$AvgSignal, data_sncrna$Functional, summary) 

# Split data by 'Functional' for t-test
negative_control <- data_sncrna[data_sncrna$Functional == "No", "AvgSignal"]
positive_data <- data_sncrna[data_sncrna$Functional == "Yes", "AvgSignal"]

# Kolmogorov-Smirnov Test (Normality)
ks.test(negative_control, "pnorm", mean(negative_control), sd(negative_control))
ks.test(positive_data, "pnorm", mean(positive_data), sd(positive_data))

# Student's t-test (assuming equal variances)
#t.test(negative_control, positive_data, var.equal = TRUE)

# Mann-Whitney U Test (aka Wilcoxon Rank-Sum Test)
wilcox.test(negative_control, positive_data)
############################################


#######################################################################
# Sample data vectors for positive and negative controls PROTEIN CODING
protein_file <- "../data/histone_feature/H3k36me3/H3k36me3_protein_matrix.csv"
data_prot <- read.csv(protein_file, header=TRUE, sep="\t")

# Calculate basic statistics
summary(data_prot$AvgSignal)
# Statistics split by 'Functional' category
by(data_prot$AvgSignal, data_prot$Functional, summary) 

# Split data by 'Functional' for t-test
negative_control <- data_prot[data_prot$Functional == "No", "AvgSignal"]
positive_data <- data_prot[data_prot$Functional == "Yes", "AvgSignal"]

# Kolmogorov-Smirnov Test (Normality)
ks.test(negative_control, "pnorm", mean(negative_control), sd(negative_control))
ks.test(positive_data, "pnorm", mean(positive_data), sd(positive_data))

# Student's t-test (assuming equal variances)
#t.test(negative_control, positive_data, var.equal = TRUE)

# Mann-Whitney U Test (aka Wilcoxon Rank-Sum Test)
wilcox.test(negative_control, positive_data)
############################################


##########################
# Fisher's Test for LNCRNA
contingency_table <- table(data_lncrna$Functional, data_lncrna$HistonePresence)
contingency_table

fisher.test(contingency_table)
##############################

##########################
# Fisher's Test for SHORT NCRNA
contingency_table <- table(data_sncrna$Functional, data_sncrna$HistonePresence)
contingency_table

fisher.test(contingency_table)
##############################

##########################
# Fisher's Test for PROTEIN CODING
contingency_table <- table(data_prot$Functional, data_prot$HistonePresence)
contingency_table

fisher.test(contingency_table)
##############################

#########################
# Define a function to calculate Spearman correlation from a dataset
calc_spearman <- function(data, indices) {
  d <- data[indices,] # Sample data based on indices
  cor.test(d[,1], d[,2], method = "spearman")$estimate
}

# Set random seed for reproducibility
set.seed(123) 

# Bootstrap with, 1000 replications
boot_prot_ci_results <- boot(data = data.frame(data_prot[,4], data_prot[,5]), 
                     statistic = calc_spearman, 
                     R = 1000)
boot_lncrna_ci_results <- boot(data = data.frame(data_lncrna[,4], data_lncrna[,5]), 
                             statistic = calc_spearman, 
                             R = 1000)
boot_shortncrna_ci_results <- boot(data = data.frame(data_sncrna[,4], data_sncrna[,5]), 
                               statistic = calc_spearman, 
                               R = 1000)
# Get confidence intervals
boot.ci(boot_prot_ci_results, type = "perc") # Use percentile-based confidence intervals 
boot.ci(boot_lncrna_ci_results, type = "perc") # Use percentile-based confidence intervals 
boot.ci(boot_shortncrna_ci_results, type = "perc") # Use percentile-based confidence intervals 



# Create a graph to visualice features
df <- data.frame(sequence_type = c("lncRNA Exons","Short ncRNAs","Protein-coding Exons"),
                 correlation = c(correlation_list[[lncrna_file1]][1,3], 
                                 correlation_list[[sncrna_file1]][1,3], 
                                 correlation_list[[protein_file1]][1,3],
                                 
                                 correlation_list[[lncrna_file2]][1,3], 
                                 correlation_list[[sncrna_file2]][1,3], 
                                 correlation_list[[protein_file2]][1,3],
                                 
                                 correlation_list[[lncrna_file3]][1,3], 
                                 correlation_list[[sncrna_file3]][1,3], 
                                 correlation_list[[protein_file3]][1,3],
                                 
                                 correlation_list[[lncrna_file4]][1,3], 
                                 correlation_list[[sncrna_file4]][1,3], 
                                 correlation_list[[protein_file4]][1,3],
                                 
                                 correlation_list[[lncrna_file5]][1,3], 
                                 correlation_list[[sncrna_file5]][1,3], 
                                 correlation_list[[protein_file5]][1,3],
                                 
                                 correlation_list[[lncrna_file6]][1,3], 
                                 correlation_list[[sncrna_file6]][1,3], 
                                 correlation_list[[protein_file6]][1,3],
                                 
                                 correlation_list[[lncrna_file7]][1,3], 
                                 correlation_list[[sncrna_file7]][1,3], 
                                 correlation_list[[protein_file7]][1,3],

                                 correlation_list[[lncrna_file8]][1,3], 
                                 correlation_list[[sncrna_file8]][1,3], 
                                 correlation_list[[protein_file8]][1,3]),
                 
                 lower_ci = c(correlation_list[[lncrna_file1]][2,3], 
                              correlation_list[[sncrna_file1]][2,3], 
                              correlation_list[[protein_file1]][2,3],
                              
                              correlation_list[[lncrna_file2]][2,3], 
                              correlation_list[[sncrna_file2]][2,3], 
                              correlation_list[[protein_file2]][2,3],
                              
                              correlation_list[[lncrna_file3]][2,3], 
                              correlation_list[[sncrna_file3]][2,3], 
                              correlation_list[[protein_file3]][2,3],
                              
                              correlation_list[[lncrna_file4]][2,3], 
                              correlation_list[[sncrna_file4]][2,3], 
                              correlation_list[[protein_file4]][2,3],
                              
                              correlation_list[[lncrna_file5]][2,3], 
                              correlation_list[[sncrna_file5]][2,3], 
                              correlation_list[[protein_file5]][2,3],
                              
                              correlation_list[[lncrna_file6]][2,3], 
                              correlation_list[[sncrna_file6]][2,3], 
                              correlation_list[[protein_file6]][2,3],
                              
                              correlation_list[[lncrna_file7]][2,3], 
                              correlation_list[[sncrna_file7]][2,3], 
                              correlation_list[[protein_file7]][2,3],

                              correlation_list[[lncrna_file8]][2,3], 
                              correlation_list[[sncrna_file8]][2,3], 
                              correlation_list[[protein_file8]][2,3]),
                 
                 upper_ci = c(correlation_list[[lncrna_file1]][3,3], 
                              correlation_list[[sncrna_file1]][3,3], 
                              correlation_list[[protein_file1]][3,3],
                              
                              correlation_list[[lncrna_file2]][3,3], 
                              correlation_list[[sncrna_file2]][3,3], 
                              correlation_list[[protein_file2]][3,3],
                              
                              correlation_list[[lncrna_file3]][3,3], 
                              correlation_list[[sncrna_file3]][3,3], 
                              correlation_list[[protein_file3]][3,3],
                              
                              correlation_list[[lncrna_file4]][3,3], 
                              correlation_list[[sncrna_file4]][3,3], 
                              correlation_list[[protein_file4]][3,3],
                              
                              correlation_list[[lncrna_file5]][3,3], 
                              correlation_list[[sncrna_file5]][3,3], 
                              correlation_list[[protein_file5]][3,3],
                              
                              correlation_list[[lncrna_file6]][3,3], 
                              correlation_list[[sncrna_file6]][3,3], 
                              correlation_list[[protein_file6]][3,3],
                              
                              correlation_list[[lncrna_file7]][3,3], 
                              correlation_list[[sncrna_file7]][3,3], 
                              correlation_list[[protein_file7]][3,3],

                              correlation_list[[lncrna_file8]][3,3], 
                              correlation_list[[sncrna_file8]][3,3], 
                              correlation_list[[protein_file8]][3,3]),
                 feature_type = c("H3k36me3", "H3k36me3", "H3k36me3",
                                  "H3k27ac", "H3k27ac", "H3k27ac",
                                  "H3k4me3", "H3k4me3", "H3k4me3",
                                  "H3k9me3", "H3k9me3", "H3k9me3",
                                  "H3k4me1", "H3k4me1", "H3k4me1",
                                  "H3k27me3", "H3k27me3", "H3k27me3",
                                  "Chromatin Accesibility","Chromatin Accesibility","Chromatin Accesibility",
                                  "Methylome","Methylome","Methylome"))

# Create a factor with the correct order for the x-axis
df$feature_name <- factor(df$sequence_type, levels = unique(df$sequence_type))

# Create a factor for sequence_type with the correct order for shapes
df$sequence_type <- factor(df$sequence_type, 
                           levels = c("lncRNA Exons","Short ncRNAs","Protein-coding Exons"))

# Create a factor for feature_type with the correct order for colors
df$feature_type <- factor(df$feature_type, 
                          levels = c("H3k36me3", "H3k27ac","H3k4me3", "H3k9me3","H3k4me1", "H3k27me3","Chromatin Accesibility","Methylome"))

ggplot(df, aes(x = feature_type, y = correlation, color = feature_type, shape = sequence_type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.05, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(19, 17, 15)) +  #c(19, 17, 15))
  scale_color_manual(values = c("blue","green","pink","orange","brown","purple","red","yellow")) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
        legend.position = "right") +
  labs(title = "Spearman correlation Histone Marks Analysis",
       x = "",
       y = "Spearman correlation",
       color = "Feature Type",
       shape = "Sequence Type")



args <- commandArgs(trailingOnly = TRUE)
#file <- args[1]


options(repos = "https://cloud.r-project.org")
install.packages("RVAideMemoire")
library(RVAideMemoire)
install.packages("randomForest")
library(randomForest)

# 1. Read the data from the file into a data frame
spearman_correlation <- function(file){
  
  data <- read.csv(file, header=TRUE, sep = "\t")
  
  # Convert Functional "Yes" "No" into numeric values 
  data$Functional <- factor(data$Functional, levels=c("No","Yes")) 
  functional_numeric <- as.numeric(factor(data$Functional)) - 1 
  functional_numeric
  # Approx NA values (CHECK method) 
  #data <- na.roughfix(data)
  
  # 2. Extract the other numeric columns
  numeric_data <- data[, !names(data) %in% c("Chromosome","ID","Functional")]
  
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


#H3k36me3
lncrna_file1 <- "../data/histone_feature/H3k36me3/H3k36me3_lncrna_matrix.csv"
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file1 <- "../data/histone_feature/H3k36me3/H3k36me3_short_ncrna_matrix.csv"
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file1 <- "../data/histone_feature/H3k36me3/H3k36me3_protein_matrix.csv"
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"


# Initialize empty list to store correlation vectors
correlation_list <- list()

correlation_list[[lncrna_file1]] <- spearman_correlation(lncrna_file1)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file1]] <- spearman_correlation(sncrna_file1)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file1]] <- spearman_correlation(protein_file1)
#correlation_matrix_combined <- do.call(rbind, correlation_list)



#H3k27ac
lncrna_file2 <- "../data/histone_feature/H3k27ac/H3k27ac_lncrna_matrix.csv"
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file2 <- "../data/histone_feature/H3k27ac/H3k27ac_short_ncrna_matrix.csv"
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file2 <- "../data/histone_feature/H3k27ac/H3k27ac_protein_matrix.csv"
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"


correlation_list[[lncrna_file2]] <- spearman_correlation(lncrna_file2)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file2]] <- spearman_correlation(sncrna_file2)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file2]] <- spearman_correlation(protein_file2)
#correlation_matrix_combined <- do.call(rbind, correlation_list)





#H3k4me3
lncrna_file3 <- "../data/histone_feature/H3k4me3/H3k4me3_lncrna_matrix.csv"
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file3 <- "../data/histone_feature/H3k4me3/H3k4me3_short_ncrna_matrix.csv"
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file3 <- "../data/histone_feature/H3k4me3/H3k4me3_protein_matrix.csv"
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"



correlation_list[[lncrna_file3]] <- spearman_correlation(lncrna_file3)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file3]] <- spearman_correlation(sncrna_file3)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file3]] <- spearman_correlation(protein_file3)
#correlation_matrix_combined <- do.call(rbind, correlation_list)




#H3k9me3
lncrna_file4 <- "../data/histone_feature/H3k9me3/H3k9me3_lncrna_matrix.csv"
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file4 <- "../data/histone_feature/H3k9me3/H3k9me3_short_ncrna_matrix.csv"
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file4 <- "../data/histone_feature/H3k9me3/H3k9me3_protein_matrix.csv"
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"


correlation_list[[lncrna_file4]] <- spearman_correlation(lncrna_file4)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file4]] <- spearman_correlation(sncrna_file4)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file4]] <- spearman_correlation(protein_file4)
#correlation_matrix_combined <- do.call(rbind, correlation_list)



#H3k4me1
lncrna_file5 <- "../data/histone_feature/H3k4me1/H3k4me1_lncrna_matrix.csv"
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file5 <- "../data/histone_feature/H3k4me1/H3k4me1_short_ncrna_matrix.csv"
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file5 <- "../data/histone_feature/H3k4me1/H3k4me1_protein_matrix.csv"
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"



correlation_list[[lncrna_file5]] <- spearman_correlation(lncrna_file5)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file5]] <- spearman_correlation(sncrna_file5)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file5]] <- spearman_correlation(protein_file5)
#correlation_matrix_combined <- do.call(rbind, correlation_list)




#H3k27me3
lncrna_file6 <- "../data/histone_feature/H3k27me3/H3k27me3_lncrna_matrix.csv"
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file6 <- "../data/histone_feature/H3k27me3/H3k27me3_short_ncrna_matrix.csv"
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file6 <- "../data/histone_feature/H3k27me3/H3k27me3_protein_matrix.csv"
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"


correlation_list[[lncrna_file6]] <- spearman_correlation(lncrna_file6)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file6]] <- spearman_correlation(sncrna_file6)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file6]] <- spearman_correlation(protein_file6)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

#chromatin accessibility
lncrna_file7 <- "../data/chrm_acc_feature/chrm_acc_lncrna_matrix.csv"
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file7 <- "../data/chrm_acc_feature/chrm_acc_short_ncrna_matrix.csv"
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file7 <- "../data/chrm_acc_feature/chrm_acc_protein_matrix.csv"
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"


correlation_list[[lncrna_file7]] <- spearman_correlation(lncrna_file7)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file7]] <- spearman_correlation(sncrna_file7)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file7]] <- spearman_correlation(protein_file7)
#correlation_matrix_combined <- do.call(rbind, correlation_list)


#methylome
lncrna_file8 <- "../data/methylome_feature/lncrna_matrix.csv"
#output_file <- "../data/histone_feature/lncrna_spearman_matrix.csv"

sncrna_file8 <- "../data/methylome_feature/short_ncrna_matrix.csv"
#output_file <- "../data/histone_feature/protein_spearman_matrix.csv"

protein_file8 <- "../data/methylome_feature/protein_matrix.csv"
#output_file <- "../data/histone_feature/short_ncrna_spearman_matrix.csv"


correlation_list[[lncrna_file8]] <- spearman_correlation(lncrna_file8)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[sncrna_file8]] <- spearman_correlation(sncrna_file8)
#correlation_matrix_combined <- do.call(rbind, correlation_list)

correlation_list[[protein_file8]] <- spearman_correlation(protein_file8)
#correlation_matrix_combined <- do.call(rbind, correlation_list)



# Write the correlation matrix to a CSV file
write.csv(correlation_matrix_combined, file = output_file, row.names = TRUE)



# Find out range for Histone marks
histone_data <- read.csv("temp_histone_marks.bed", header=FALSE, sep = "\t")
summary(histone_data$V7)

histone_data <- read.csv("../data/histone_feature/lncrna-histone-feature-matrix.csv", header=TRUE)
summary(histone_data$AvgSignal)


#######################################
# Statistics for methylation files:
methyl_file <- "../data/ENCFF913ZNZ.bed"
methyl_data <- read.csv(methyl_file, header = FALSE, sep = "\t", nrows = 100000)
summary(methyl_data)

