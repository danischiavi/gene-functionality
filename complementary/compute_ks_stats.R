# Generate a heatmap of gene functionality feature correlation.

# Spearman Correlation Heatmap plot
# Load necessary libraries
library(corrr)
library(Hmisc) # For rcorr() function
library(ggcorrplot)
library(dplyr)
library(reshape2)

# Load KS statistics utility functions. (A custom version of ks to keep the sign in the D statistic)
source("complementary/config.R")
source("complementary/utils.R")

# Load data
source("complementary/load_gene_functionality_features.R")
feature_matrix <- load_gene_functionality_features()

# Set features to analyze
HEATMAP_SELECT_FEATURES <- c("Random number","GC content", "low_complexity_density", 
                             "CpG", "GA", "TA", "GT",
                             "phyloP max_241w", "phyloP max_100w", 
                             "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                             "RPKM_tissue", "RPKM_primary cell", 
                             "H3K9ac", "H3K79me1", "H3K79me2", 
                             "chromatin_acc", "methylome", 
                             "Repeat free", "Copy number", "RNAcode", 
                             "Fickett_score", "Max covariance", "MFE",
                             "Interaction_ave", "gnomAD_SNP_density", "gnomAD_MAF",
                             "T content", "C content", "G content", "MaxEnt splicing")

# Set features to analyze
colnames(feature_matrix)
HEATMAP_SELECT_FEATURES_OLD_NAMES <- c("Random","GC_percentage", "lowComplexity_density", 
                                       "CpG", "GA", "GT", "TA",
                                       "phyloP_max_241w", "phyloP_max_100w", 
                                       "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                                       "RPKM_tissue", "RPKM_primary.cell", 
                                       "H3K9ac_MaxScaledSignal", "H3K79me1_MaxScaledSignal", "H3K79me2_MaxScaledSignal", 
                                       "chrm_acc_MaxScaledSignal", "methylome", 
                                       "repeat_distance", "copy_number", "coding_potential", 
                                       "fickett", "Max_covariance", "MFE",
                                       "Interaction_ave", "SNP_density", "MAF_avg",
                                       "T", "C", "G", "MaxEnt_splicing")

HEATMAP_SELECT_FEATURES_LABELS <- c("Random number", "GC%", "Complexity",
                                    "CpG",  "GA", "GT", "TA",
                                    "PhyloP-mammals", "PhyloP-vertebrates",
                                    "GERP-mammals", "GERP-vertebrates",
                                    "Tissue RPKM", "Primary cell RPKM", 
                                    "H3K9ac", "H3K79me1", "H3K79me2",
                                    "Chromatin Accessibility", "Methylome",
                                    "Repeat free", "Copies", "RNAcode",
                                    "Fickett", "Covariance", "MFE",
                                    "Interactions", "SNPs", "MAF",
                                    "T%", "C%", "G%", "MaxEnt splicing")

# Select features and convert to numeric
data_numeric <- feature_matrix %>% 
  dplyr::select(all_of(HEATMAP_SELECT_FEATURES_OLD_NAMES)) %>%
  sapply(function(feature) as.numeric(as.character(feature))) %>%
  as.data.frame()
# Restore Dataset column
data_numeric$Dataset <- feature_matrix$Dataset


######################################
# Compute ks-stats for raw features. #
########
# mrna #
subset_m <- data_numeric %>%
  filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>%
  dplyr::select(-Dataset)
subset_nc_m <- data_numeric %>%
  filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>%
  dplyr::select(-Dataset)

# Run the K-S tests
ks_results_m <- run_ks_tests(subset_nc_m, subset_m)

ks_d_values_m <- as.data.frame(t(ks_results_m))
rownames(ks_d_values_m) <- HEATMAP_SELECT_FEATURES_LABELS

ks_pvalues_m <- ks_d_values_m[["p.val"]]
# Adjust p-values using Benjamini
ks_pvalues_adj_m <- p.adjust(ks_pvalues_m, method="BH")
ks_d_values_m$p.val.adj = ks_pvalues_adj_m

# Save results
write.csv(ks_d_values_m, GENE_FUNCTIONALITY_MRNA_KS_STAT_FILE)

##########
# lncrna #
subset_l <- data_numeric %>%
  filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")) %>%
  dplyr::select(-Dataset)
subset_nc_l <- data_numeric %>%
  filter(Dataset %in% c("lncrna-exon1-negative-control","lncrna-exon2-negative-control")) %>%
  dplyr::select(-Dataset)

# Run the K-S tests
ks_results_l <- run_ks_tests(subset_nc_l, subset_l)

ks_d_values_l <- as.data.frame(t(ks_results_l))
rownames(ks_d_values_l) <- HEATMAP_SELECT_FEATURES_LABELS

ks_pvalues_l <- ks_d_values_l[["p.val"]]
# Adjust p-values using Benjamini
ks_pvalues_adj_l <- p.adjust(ks_pvalues_l, method="BH")
ks_d_values_l$p.val.adj = ks_pvalues_adj_l

#print(ks_results_l)
write.csv(ks_d_values_l, GENE_FUNCTIONALITY_LNCRNA_KS_STAT_FILE)

##########
# sncrna #
subset_s <- data_numeric %>%
  filter(Dataset %in% c("short-ncrna")) %>%
  dplyr::select(-Dataset)
subset_nc_s <- data_numeric %>%
  filter(Dataset %in% c("short-ncrna-negative-control")) %>%
  dplyr::select(-Dataset)

# Run the K-S tests
ks_results_s <- run_ks_tests(subset_nc_s, subset_s)

ks_d_values_s <- as.data.frame(t(ks_results_s))
rownames(ks_d_values_s) <- HEATMAP_SELECT_FEATURES_LABELS

ks_pvalues_s <- ks_d_values_s[["p.val"]]
# Adjust p-values using Benjamini
ks_pvalues_adj_s <- p.adjust(ks_pvalues_s, method="BH")
ks_d_values_s$p.val.adj = ks_pvalues_adj_s

#print(ks_results_s)
write.csv(ks_d_values_s, GENE_FUNCTIONALITY_SNCRNA_KS_STAT_FILE)
