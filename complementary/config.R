# Configuration file for project paths

# Base directory for data
DATA_DIR <- "data"

# Complementary directory
COMPLEMENTARY_DIR <- "complementary"

# Base directory for results
RESULTS_DIR <- "results"

# Subdirectories within results
ZSCORES_DIR <- file.path(RESULTS_DIR, "z-scores")

KS_STATS_DIR <- file.path(COMPLEMENTARY_DIR, "ks_stats")

DISTANCE_EFFECT_DIR <- file.path(RESULTS_DIR, "Paper figures")
DISTANCE_EFFECT_HEATSCATTER_DIR <- DISTANCE_EFFECT_DIR
DISTANCE_EFFECT_SPEARMAN_DIR <- DISTANCE_EFFECT_DIR

VIOLIN_PLOTS_DIR <- file.path(RESULTS_DIR, "Paper figures")
DISTRIBUTION_PLOTS_DIR <- file.path(VIOLIN_PLOTS_DIR, "Distribution_plots")

LATEST_LOG_SCATTER_PLOTS_DIR <- file.path(RESULTS_DIR, "scatter_plots")


# Specific file paths
# Data files
EPIGENETIC_HISTONES_DIR <- file.path(COMPLEMENTARY_DIR, "histone_marks")

#FEATURES_20260608_DIR <- file.path(FEATURES_DIR, "20260608")

FUNC_PROT_EXON2_FEATURES_FILE <- file.path(RESULTS_DIR, "functional-protein-exon2-dataset-features_20260608.csv")
FUNC_PROT_EXON3_FEATURES_FILE <- file.path(RESULTS_DIR, "functional-protein-exon3-dataset-features_20260608.csv")
FUNC_LNCRNA_EXON1_FEATURES_FILE <- file.path(RESULTS_DIR, "functional-lncrna-exon1-dataset-features_20260608.csv")
FUNC_LNCRNA_EXON2_FEATURES_FILE <- file.path(RESULTS_DIR, "functional-lncrna-exon2-dataset-features_20260608.csv")
FUNC_SNCRNA_FEATURES_FILE <- file.path(RESULTS_DIR, "functional-short-ncrna-dataset-features_20260608.csv")

NC_PROT_EXON2_FEATURES_FILE <- file.path(RESULTS_DIR, "protein-exon2-negative-control-dataset-features_20260608.csv")
NC_PROT_EXON3_FEATURES_FILE <- file.path(RESULTS_DIR, "protein-exon3-negative-control-dataset-features_20260608.csv")
NC_LNCRNA_EXON1_FEATURES_FILE <- file.path(RESULTS_DIR, "lncrna-exon1-negative-control-dataset-features_20260608.csv")
NC_LNCRNA_EXON2_FEATURES_FILE <- file.path(RESULTS_DIR, "lncrna-exon2-negative-control-dataset-features_20260608.csv")
NC_SNCRNA_FEATURES_FILE <- file.path(RESULTS_DIR, "short-ncrna-negative-control-dataset-features_20260608.csv")


# Epigenetic histone files
COMBINED_PROT_EXON2_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_protein-exon2.csv")
COMBINED_PROT_EXON3_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_protein-exon3.csv")
COMBINED_LNCRNA_EXON1_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_lncrna-exon1.csv")
COMBINED_LNCRNA_EXON2_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_lncrna-exon2.csv")
COMBINED_SNCRNA_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_short-ncrna.csv")

COMBINED_NC_PROT_EXON2_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_protein-exon2-NC.csv")
COMBINED_NC_PROT_EXON3_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_protein-exon3-NC.csv")
COMBINED_NC_LNCRNA_EXON1_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_lncrna-exon1-NC.csv")
COMBINED_NC_LNCRNA_EXON2_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_lncrna-exon2-NC.csv")
COMBINED_NC_SNCRNA_HISTONE_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "combined_short-ncrna-NC.csv")

# Epigenetic chromatin files
CHR_ACC_PROT_EXON2_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_protein-exon2-chrm_acc-feature.csv")
CHR_ACC_PROT_EXON3_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_protein-exon3-chrm_acc-feature.csv")
CHR_ACC_LNCRNA_EXON1_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_lncrna-exon1-chrm_acc-feature.csv")
CHR_ACC_LNCRNA_EXON2_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_lncrna-exon2-chrm_acc-feature.csv")
CHR_ACC_SNCRNA_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_short-ncrna-chrm_acc-feature.csv")

NC_CHR_ACC_PROT_EXON2_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_protein-exon2-NC-chrm_acc-feature.csv")
NC_CHR_ACC_PROT_EXON3_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_protein-exon3-NC-chrm_acc-feature.csv")
NC_CHR_ACC_LNCRNA_EXON1_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_lncrna-exon1-NC-chrm_acc-feature.csv")
NC_CHR_ACC_LNCRNA_EXON2_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_lncrna-exon2-NC-chrm_acc-feature.csv")
NC_CHR_ACC_SNCRNA_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "chrm_acc_short-ncrna-NC-chrm_acc-feature.csv")


# Methylome files
METHYLOME_PROT_EXON2_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "protein-exon2-methylome-feature.csv")
METHYLOME_PROT_EXON3_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "protein-exon3-methylome-feature.csv")
METHYLOME_LNCRNA_EXON1_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "lncrna-exon1-methylome-feature.csv")
METHYLOME_LNCRNA_EXON2_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "lncrna-exon2-methylome-feature.csv")
METHYLOME_SNCRNA_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "short-ncrna-methylome-feature.csv")

NC_METHYLOME_PROT_EXON2_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "protein-exon2-NC-methylome-feature.csv")
NC_METHYLOME_PROT_EXON3_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "protein-exon3-NC-methylome-feature.csv")
NC_METHYLOME_LNCRNA_EXON1_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "lncrna-exon1-NC-methylome-feature.csv")
NC_METHYLOME_LNCRNA_EXON2_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "lncrna-exon2-NC-methylome-feature.csv")
NC_METHYLOME_SNCRNA_FILE <- file.path(EPIGENETIC_HISTONES_DIR, "short-ncrna-NC-methylome-feature.csv")

# DAF classification file
DAF_DIR <- file.path(COMPLEMENTARY_DIR, "daf")
DAF_SNCRNA_POS_FILE <- file.path(DAF_DIR, "sncrna_positives.tsv")
DAF_SNCRNA_CON_FILE <- file.path(DAF_DIR, "sncrna_controls.tsv")


# KS Stats files
EPIGENETIC_MRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "epigenetic_mrna_ks_stat_sample.csv")
EPIGENETIC_LNCRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "epigenetic_lncrna_ks_stat_sample.csv")
EPIGENETIC_SNCRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "epigenetic_sncrna_ks_stat_sample.csv")
EPIGENETIC_MEAN_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "epigenetic_mean_ks_stat_sample.csv")
INTRINSIC_MRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "intrinsic_mrna_ks_stat_sample.csv")
INTRINSIC_LNCRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "intrinsic_lncrna_ks_stat_sample.csv")
INTRINSIC_SNCRNA_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "intrinsic_sncrna_ks_stat_sample.csv")
INTRINSIC_MEAN_KS_STAT_SAMPLE_FILE <- file.path(KS_STATS_DIR, "intrinsic_mean_ks_stat_sample.csv")
GENE_FUNCTIONALITY_MRNA_KS_STAT_FILE <- file.path(KS_STATS_DIR, "gene_functionality_mrna_ks_stat_sample.csv")
GENE_FUNCTIONALITY_LNCRNA_KS_STAT_FILE <- file.path(KS_STATS_DIR, "gene_functionality_lncrna_ks_stat_sample.csv")
GENE_FUNCTIONALITY_SNCRNA_KS_STAT_FILE <- file.path(KS_STATS_DIR, "gene_functionality_sncrna_ks_stat_sample.csv")
GENE_FUNCTIONALITY_MEAN_KS_STAT_FILE <- file.path(KS_STATS_DIR, "gene_functionality_mean_ks_stat_sample.csv")

# Distance effect files
EFFECT_ON_DISTANCE_JOINED_FILE <- file.path(DISTANCE_EFFECT_HEATSCATTER_DIR, "effectOnDistanceHeatscatterJoined5k_5M.png")
CORR_MATRIX_PROTEIN_FILE <- file.path(DISTANCE_EFFECT_SPEARMAN_DIR, "corr_matrix_protein.csv")
CORR_MATRIX_SNCRNA_FILE <- file.path(DISTANCE_EFFECT_SPEARMAN_DIR, "corr_matrix_sncrna.csv")
CORR_MATRIX_LNCRNA_FILE <- file.path(DISTANCE_EFFECT_SPEARMAN_DIR, "corr_matrix_lncrna.csv")


# Violin plot files (legacy subset paths — VIOLIN_PLOTS_SUBSET_DIR not yet defined)
#EPIGENETIC_P1_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p1_1.png")
#EPIGENETIC_P1_2_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p1_2.png")
#EPIGENETIC_P2_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p2_1.png")
#EPIGENETIC_P2_2_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p2_2.png")
#EPIGENETIC_P3_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p3_1.png")
#EPIGENETIC_P4_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p4_1.png")
#EPIGENETIC_P5_1_PLOT_FILE <- file.path(VIOLIN_PLOTS_SUBSET_DIR, "Epigenetic_p5_1.png")


# PCA plot files (legacy — PCA_DIR not yet defined)
#PCA_PROTEIN_20_FEATURES_PLOT_FILE <- file.path(PCA_DIR, "pca_protein_20_features.png")
#PCA_LNCRNA_20_FEATURES_PLOT_FILE <- file.path(PCA_DIR, "pca_lncrna_20_features.png")
#PCA_SNCRNA_20_FEATURES_PLOT_FILE <- file.path(PCA_DIR, "pca_sncrna_20_features.png")
#PCA_20_FEATURES_JOINED_PLOT_FILE <- file.path(PCA_DIR, "pca_20_features_joined.png")
#PCA_PROTEIN_20_FEATURES_LOADINGS_PLOT_FILE <- file.path(PCA_DIR, "pca_protein_20_features_loadings.png")
#PCA_LNCRNA_20_FEATURES_LOADINGS_PLOT_FILE <- file.path(PCA_DIR, "pca_lncrna_20_features_loadings.png")
#PCA_SNCRNA_20_FEATURES_LOADINGS_PLOT_FILE <- file.path(PCA_DIR, "pca_sncrna_20_features_loadings.png")
#PCA_20_FEATURES_LOADINGS_JOINED_PLOT_FILE <- file.path(PCA_DIR, "pca_20_features_loadings_joined.png")

