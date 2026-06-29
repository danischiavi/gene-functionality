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

dir.create(file.path(VIOLIN_PLOTS_DIR, "png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(VIOLIN_PLOTS_DIR, "pdf"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(KS_STATS_DIR), showWarnings = FALSE, recursive = TRUE)
dir.create("complementary/spearman_correlations/gene_features", showWarnings = FALSE, recursive = TRUE)
dir.create("complementary/spearman_correlation", showWarnings = FALSE, recursive = TRUE)

# Load data
source("complementary/load_gene_functionality_features.R")
feature_matrix <- load_gene_functionality_features()

# Set features to analyze
HEATMAP_SELECT_FEATURES <- c("Random number","GC content", "low_complexity_density", 
                             "CpG", "GA", "TA", "GT",
                             "T", "C", "G",
                             "phyloP max_241w", "phyloP max_100w", 
                             "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                             "RPKM_tissue", "RPKM_primary cell", 
                             "H3K9ac", "H3K79me1", "H3K79me2", 
                             "chromatin_acc", "methylome", 
                             "Repeat free", "Copy number", "RNAcode", 
                             "Fickett_score", "Max covariance", "MFE",
                             "Interaction_ave", "MaxEnt_splicing", "gnomAD_SNP_density", "gnomAD_MAF")

# Set features to analyze
colnames(feature_matrix)
HEATMAP_SELECT_FEATURES_OLD_NAMES <- c("Random","GC_percentage", "lowComplexity_density", 
                             "CpG", "GA", "GT", "TA",
                             "T", "C", "G",
                             "phyloP_max_241w", "phyloP_max_100w", 
                             "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                             "RPKM_tissue", "RPKM_primary.cell", 
                             "H3K9ac_MaxScaledSignal", "H3K79me1_MaxScaledSignal", "H3K79me2_MaxScaledSignal", 
                             "chrm_acc_MaxScaledSignal", "methylome", 
                             "repeat_distance", "copy_number", "coding_potential", 
                             "fickett", "Max_covariance", "MFE",
                             "Interaction_ave", "MaxEnt_splicing", "SNP_density", "MAF_avg")

HEATMAP_SELECT_FEATURES_LABELS <- c("Random number", "GC%", "Complexity",
                                    "CpG",  "GA", "GT", "TA",
                                    "T%", "C%", "G%",
                                    "PhyloP-mammals", "PhyloP-vertebrates",
                                    "GERP-mammals", "GERP-vertebrates",
                                    "Tissue RPKM", "Primary cell RPKM", 
                                    "H3K9ac", "H3K79me1", "H3K79me2",
                                    "Chromatin Accessibility", "Methylome",
                                    "Repeat free", "Copies", "RNAcode",
                                    "Fickett", "Covariance", "MFE",
                                    "Interactions", "MaxEnt splicing", "SNPs", "MAF")

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

##################
# mean|signed_D| #
ks_results <- data.frame(
  meanAbsD = (abs(ks_d_values_m$signed_D) + abs(ks_d_values_s$signed_D) + abs(ks_d_values_l$signed_D))/3
)
rownames(ks_results) <- HEATMAP_SELECT_FEATURES_LABELS
#print(ks_results)
write.csv(ks_results, GENE_FUNCTIONALITY_MEAN_KS_STAT_FILE)


###############################
# Compute Spearman Correlations
prot_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3",
                            "protein-exon2-negative-control", "protein-exon3-negative-control")) %>%
      dplyr::select(all_of(HEATMAP_SELECT_FEATURES_OLD_NAMES))
      #select(all_of(SELECTED_HISTONE_FEATURES))
    ),
  type = "spearman"
)


lncrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2",
                            "lncrna-exon1-negative-control", "lncrna-exon2-negative-control")) %>%
      dplyr::select(all_of(HEATMAP_SELECT_FEATURES_OLD_NAMES))
      #select(all_of(SELECTED_HISTONE_FEATURES))
  ),
  type = "spearman"
)

sncrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset %in% c("short-ncrna", "short-ncrna-negative-control")) %>%
      dplyr::select(all_of(HEATMAP_SELECT_FEATURES_OLD_NAMES))
      #dplyr::select(all_of(SELECTED_HISTONE_FEATURES))
  ),
  type = "spearman"
)

# Protein coding #
prot_corr_matrix <- prot_corr_obj[["r"]]
prot_corr_matrix[is.na(prot_corr_matrix)] <- NA
prot_pvalues <- prot_corr_obj[["P"]]

# Number of features
n <- ncol(prot_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Benjamini-H
prot_p_values_adjusted <- p.adjust(prot_pvalues, method="BH")

prot_p_values_adjusted <- matrix(prot_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(prot_p_values_adjusted) <- diag

rownames(prot_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(prot_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS

rownames(prot_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(prot_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS


# Generate heatmap.
ggcorrplot(prot_corr_matrix, type = "lower", method = "square", lab = TRUE, 
           lab_col = "black", lab_size = 2.5, ggtheme = theme_void, 
           title = "Spearman correlation Heatmap - mRNA", show.diag = TRUE, 
           p.mat = prot_p_values_adjusted, sig.level = corrected_alpha, insig = "pch"
) +
  theme(
    plot.title = element_text(size = 33),  # Increase title size
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    
  )

# Save heatmap 
ggsave("figureS4-mRNA-correlation-matrix.png", path = file.path(VIOLIN_PLOTS_DIR, "png"), scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
ggsave("figureS4-mRNA-correlation-matrix.pdf", path = file.path(VIOLIN_PLOTS_DIR, "pdf"), scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

write.csv(prot_corr_matrix, "complementary/spearman_correlations/gene_features/protein-correlation-matrix_p&n_20260611.csv",)
write.csv(prot_p_values_adjusted, "complementary/spearman_correlations/gene_features/protein-correlation-p-adj_p&n_20260611.csv")

protein_corr_matrix_melted <- melt(prot_corr_matrix)
protein_p_values_adjusted_melted <- melt(prot_p_values_adjusted)

write.csv(protein_corr_matrix_melted, "complementary/spearman_correlations/gene_features/protein-correlation-matrix_melted_p&n_20260611.csv")
write.csv(protein_p_values_adjusted_melted, "complementary/spearman_correlations/gene_features/protein-correlation-p-adj_melt_p&n_20260611.csv")


# lncRNA #
lncrna_corr_matrix <- lncrna_corr_obj[["r"]]
lncrna_pvalues <- lncrna_corr_obj[["P"]]

# Number of features
n <- ncol(lncrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
lncrna_p_values_adjusted <- p.adjust(lncrna_pvalues, method="BH")

lncrna_p_values_adjusted <- matrix(lncrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(lncrna_p_values_adjusted) <- diag

rownames(lncrna_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(lncrna_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
#rownames(lncrna_corr_matrix) <- SELECTED_HISTONE_LABELS
#colnames(lncrna_corr_matrix) <- SELECTED_HISTONE_LABELS

rownames(lncrna_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(lncrna_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
#rownames(lncrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS
#colnames(lncrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS

ggcorrplot(lncrna_corr_matrix, type = "lower", method = "square", lab = TRUE, 
           lab_col = "black", lab_size = 2.5, ggtheme = theme_void, 
           title = "Spearman correlation Heatmap - lncRNA", show.diag = TRUE, 
           p.mat = lncrna_p_values_adjusted, sig.level = corrected_alpha, insig = "pch"
) +
  theme(
    plot.title = element_text(size = 33),  # Increase title size
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    
  )
ggsave("figureS6-lncRNA-correlation-matrix.png", path = file.path(VIOLIN_PLOTS_DIR, "png"), scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
ggsave("figureS6-lncRNA-correlation-matrix.pdf", path = file.path(VIOLIN_PLOTS_DIR, "pdf"), scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

write.csv(lncrna_corr_matrix, "complementary/spearman_correlations/gene_features/lncrna-correlation-matrix_p&n_20260611.csv")
write.csv(lncrna_p_values_adjusted, "complementary/spearman_correlations/gene_features/lncrna-correlation-p-adj_p&n_20260611.csv")

lncrna_corr_matrix_melted <- melt(lncrna_corr_matrix)
lncrna_p_values_adjusted_melted <- melt(lncrna_p_values_adjusted)

write.csv(lncrna_corr_matrix_melted, "complementary/spearman_correlations/gene_features/lncrna-correlation-matrix_melted_p&n_20260611.csv")
write.csv(lncrna_p_values_adjusted_melted, "complementary/spearman_correlations/gene_features/lncrna-correlation-p-adj_melt_p&n_20260611.csv")


# sncRNA #
sncrna_corr_matrix <- sncrna_corr_obj[["r"]]
sncrna_pvalues <- sncrna_corr_obj[["P"]]

# Number of features
n <- ncol(sncrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
sncrna_p_values_adjusted <- p.adjust(sncrna_pvalues, method="BH")

sncrna_p_values_adjusted <- matrix(sncrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(sncrna_p_values_adjusted) <- diag

rownames(sncrna_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(sncrna_corr_matrix) <- HEATMAP_SELECT_FEATURES_LABELS

rownames(sncrna_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS
colnames(sncrna_p_values_adjusted) <- HEATMAP_SELECT_FEATURES_LABELS

# Generate heatmap.
ggcorrplot(sncrna_corr_matrix, type = "lower", method = "square", lab = TRUE, 
           lab_col = "black", lab_size = 2.5, ggtheme = theme_void, 
           title = "Spearman correlation Heatmap - sncRNA", show.diag = TRUE, 
           p.mat = sncrna_p_values_adjusted, sig.level = corrected_alpha, insig = "pch"
) +
  theme(
    plot.title = element_text(size = 33),  # Increase title size
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    
  )
ggsave("figureS5-sncRNA-correlation-matrix.png", path = file.path(VIOLIN_PLOTS_DIR, "png"), scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
ggsave("figureS5-sncRNA-correlation-matrix.pdf", path = file.path(VIOLIN_PLOTS_DIR, "pdf"), scale = 3,  width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

write.csv(sncrna_corr_matrix, "complementary/spearman_correlations/gene_features/sncrna-correlation-matrix_p&n_20260611.csv")
write.csv(sncrna_p_values_adjusted, "complementary/spearman_correlations/gene_features/sncrna-correlation-p-adj_p&n_20260611.csv")

sncrna_corr_matrix_melted <- melt(sncrna_corr_matrix)
sncrna_p_values_adjusted_melted <- melt(sncrna_p_values_adjusted)

write.csv(sncrna_corr_matrix_melted, "complementary/spearman_correlations/gene_features/sncrna-correlation-matrix_melted_p&n_20260611.csv")
write.csv(sncrna_p_values_adjusted_melted, "complementary/spearman_correlations/gene_features/sncrna-correlation-p-adj_melt_p&n_20260611.csv")




########################################################
#Obtain correlation between all features just for positive cases.
# Spearman Correlation Heat map
unique(feature_matrix$Dataset)
prot_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset == "protein-coding-exon2" | Dataset == "protein-coding-exon3") %>%
      select(-Dataset)
  ),
  type = "spearman"
)


lncrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset == "lncrna-exon1" | Dataset == "lncrna-exon2") %>%
      select(-Dataset)
  ),
  type = "spearman"
)

sncrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric %>%
      filter(Dataset == "short-ncrna") %>%
      select(-Dataset)
  ),
  type = "spearman"
)

# Protein coding #
prot_corr_matrix <- prot_corr_obj[["r"]]
prot_pvalues <- prot_corr_obj[["P"]]

# Number of features
n <- ncol(prot_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
prot_p_values_adjusted <- p.adjust(prot_pvalues, method="BH")

prot_p_values_adjusted <- matrix(prot_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(prot_p_values_adjusted) <- diag

# Melt
protein_corr_matrix_melted <- melt(prot_corr_matrix)
protein_p_values_adjusted_melted <- melt(prot_p_values_adjusted)

write.csv(protein_corr_matrix_melted, "complementary/spearman_correlation/protein-correlation-matrix_melted.csv")
write.csv(protein_p_values_adjusted_melted, "complementary/spearman_correlation/protein-correlation-p-adj_melt.csv")


# lncRNA #
lncrna_corr_matrix <- lncrna_corr_obj[["r"]]
lncrna_pvalues <- lncrna_corr_obj[["P"]]

# Number of features
n <- ncol(lncrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
lncrna_p_values_adjusted <- p.adjust(lncrna_pvalues, method="BH")

lncrna_p_values_adjusted <- matrix(lncrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(lncrna_p_values_adjusted) <- diag

#Melt
lncrna_corr_matrix_melted <- melt(lncrna_corr_matrix)
lncrna_p_values_adjusted_melted <- melt(lncrna_p_values_adjusted)

write.csv(lncrna_corr_matrix_melted, "complementary/spearman_correlation/lncrna-correlation-matrix_melted.csv")
write.csv(lncrna_p_values_adjusted_melted, "complementary/spearman_correlation/lncrna-correlation-p-adj_melt.csv")


# sncRNA #
sncrna_corr_matrix <- sncrna_corr_obj[["r"]]
sncrna_pvalues <- sncrna_corr_obj[["P"]]

# Number of features
n <- ncol(sncrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Bonferroni
sncrna_p_values_adjusted <- p.adjust(sncrna_pvalues, method="BH")

sncrna_p_values_adjusted <- matrix(sncrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(sncrna_p_values_adjusted) <- diag

sncrna_corr_matrix_melted <- melt(sncrna_corr_matrix)
sncrna_p_values_adjusted_melted <- melt(sncrna_p_values_adjusted)

write.csv(sncrna_corr_matrix_melted, "complementary/spearman_correlation/sncrna-correlation-matrix_melted.csv")
write.csv(sncrna_p_values_adjusted_melted, "complementary/spearman_correlation/sncrna-correlation-p-adj_melt.csv")

