# Generate a heatmap of epigenetic feature correlation using a subset of epigenetic
# data. This script selects 100 positive cases and 1000 negative controls for each
# RNA type and their respective exons analyzed. To replicate the plots published in 
# the paper, all epigenetic data should be selected (i.e. no sub-setting).

# Load necessary libraries
library(Hmisc) # For rcorr() function
library(ggcorrplot)
library(dplyr)
library(reshape2)

library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
#install.packages("ggplotify")
library(ggplotify) ## to convert pheatmap to ggplot2
#install.packages("heatmaply")
library(heatmaply) ## for constructing interactive heatmap
library(scales)
library(plotly)       # For combining the final plots
library(htmlwidgets)  # For saving plotly as HTML
library(webshot2)     # For converting HTML to PNG/PDF (requires Chrome)


#######################################
# Load KS statistics utility functions. (A custom version of ks to keep the sign in the D statistic)
source("complementary/config.R")
source("complementary/utils.R")

source("complementary/load_epigenetic_features.R")
feature_matrix_epigenetic <- load_epigenetic_features()

###############################
# Define features of interest #
SELECTED_HISTONE_FEATURES <- c("H2AFZ_MaxScaledSignal", "H2AK5ac_MaxScaledSignal", "H2AK9ac_MaxScaledSignal",
                               "H2BK120ac_MaxScaledSignal", "H2BK12ac_MaxScaledSignal", "H2BK15ac_MaxScaledSignal",
                               "H2BK20ac_MaxScaledSignal", "H2BK5ac_MaxScaledSignal", "H3F3A_MaxScaledSignal", 
                               "H3K14ac_MaxScaledSignal",
                               "H3K18ac_MaxScaledSignal", "H3K23ac_MaxScaledSignal", "H3K23me2_MaxScaledSignal",
                               "H3K27ac_MaxScaledSignal", "H3K27me3_MaxScaledSignal", "H3K36me3_MaxScaledSignal",
                               "H3K4ac_MaxScaledSignal", "H3K4me1_MaxScaledSignal", "H3K4me2_MaxScaledSignal",
                               "H3K4me3_MaxScaledSignal", "H3K56ac_MaxScaledSignal", "H3K79me1_MaxScaledSignal",
                               "H3K79me2_MaxScaledSignal", "H3K9ac_MaxScaledSignal", "H3K9me1_MaxScaledSignal",
                               "H3K9me2_MaxScaledSignal", "H3K9me3_MaxScaledSignal", #"H3T11ph_MaxScaledSignal",
                               "H4K12ac_MaxScaledSignal", "H4K20me1_MaxScaledSignal", "H4K5ac_MaxScaledSignal",
                               "H4K8ac_MaxScaledSignal", "H4K91ac_MaxScaledSignal", "chrm_acc_MaxScaledSignal",
                               "methylome")

SELECTED_HISTONE_LABELS <- c("H2AFZ", "H2AK5ac", "H2AK9ac",
                             "H2BK120ac", "H2BK12ac", "H2BK15ac",
                             "H2BK20ac", "H2BK5ac", "H3F3A",
                             "H3K14ac",
                             "H3K18ac", "H3K23ac", "H3K23me2",
                             "H3K27ac", "H3K27me3", "H3K36me3",
                             "H3K4ac", "H3K4me1", "H3K4me2",
                             "H3K4me3", "H3K56ac", "H3K79me1",
                             "H3K79me2", "H3K9ac", "H3K9me1",
                             "H3K9me2", "H3K9me3", #"H3T11ph",
                             "H4K12ac", "H4K20me1", "H4K5ac",
                             "H4K8ac", "H4K91ac", "chrm_acc",
                             "methylome")

# Select epigenetic features and convert to numeric
data_numeric_epigenetic <- feature_matrix_epigenetic %>% 
  dplyr::select(all_of(SELECTED_HISTONE_FEATURES)) %>%
  sapply(function(feature) as.numeric(as.character(feature))) %>%
  as.data.frame()
# Restore Dataset column
data_numeric_epigenetic$Dataset <- feature_matrix_epigenetic$Dataset

# Select a subset per RNA type. Choose 100 for positive cases, 1000 for negative controls.
set.seed(2025)  # Set seed for reproducibility
data_numeric_epigenetic_sample <- rbind(
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("protein-coding-exon2")) %>%
    sample_n(100),
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("protein-coding-exon3")) %>%
    sample_n(100),
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("lncrna-exon1")) %>%
    sample_n(100),
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("lncrna-exon2")) %>%
    sample_n(100),
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("short-ncrna")) %>%
    sample_n(100),
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("protein-exon2-negative-control")) %>%
    sample_n(1000),
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("protein-exon3-negative-control")) %>%
    sample_n(1000),
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("lncrna-exon1-negative-control")) %>%
    sample_n(1000),
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("lncrna-exon2-negative-control")) %>%
    sample_n(1000),
  data_numeric_epigenetic %>%
    filter(Dataset %in% c("short-ncrna-negative-control")) %>%
    sample_n(1000)
)
# If the paper plot is desired, un-comment the following line, no other change needed.
#data_numeric_epigenetic_sample <- data_numeric_epigenetic

######################################
# Compute ks-stats for raw features. #
########
# mrna #
subset_m <- data_numeric_epigenetic_sample %>%
  filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>%
  dplyr::select(-Dataset)
subset_nc_m <- data_numeric_epigenetic_sample %>%
  filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>%
  dplyr::select(-Dataset)

# Run the K-S tests
ks_results_m <- run_ks_tests(subset_nc_m, subset_m)

#print(ks_results_m)
write.csv(ks_results_m, EPIGENETIC_MRNA_KS_STAT_SAMPLE_FILE)

ks_d_values_m <- as.data.frame(t(ks_results_m))
rownames(ks_d_values_m) <- SELECTED_HISTONE_LABELS

##########
# lncrna #
subset_l <- data_numeric_epigenetic_sample %>%
  filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")) %>%
  dplyr::select(-Dataset)
subset_nc_l <- data_numeric_epigenetic_sample %>%
  filter(Dataset %in% c("lncrna-exon1-negative-control","lncrna-exon2-negative-control")) %>%
  dplyr::select(-Dataset)

# Run the K-S tests
ks_results_l <- run_ks_tests(subset_nc_l, subset_l)

#print(ks_results_l)
write.csv(ks_results_l, EPIGENETIC_LNCRNA_KS_STAT_SAMPLE_FILE)

ks_d_values_l <- as.data.frame(t(ks_results_l))
rownames(ks_d_values_l) <- SELECTED_HISTONE_LABELS

##########
# sncrna #
subset_s <- data_numeric_epigenetic_sample %>%
  filter(Dataset %in% c("short-ncrna")) %>%
  dplyr::select(-Dataset)
subset_nc_s <- data_numeric_epigenetic_sample %>%
  filter(Dataset %in% c("short-ncrna-negative-control")) %>%
  dplyr::select(-Dataset)

# Run the K-S tests
ks_results_s <- run_ks_tests(subset_nc_s, subset_s)

#print(ks_results_s)
write.csv(ks_results_s, EPIGENETIC_SNCRNA_KS_STAT_SAMPLE_FILE)

ks_d_values_s <- as.data.frame(t(ks_results_s))
rownames(ks_d_values_s) <- SELECTED_HISTONE_LABELS

##################
# mean|signed_D| #
ks_results <- data.frame(
  meanAbsD = (abs(ks_d_values_m$signed_D) + abs(ks_d_values_s$signed_D) + abs(ks_d_values_l$signed_D))/3
)
rownames(ks_results) <- SELECTED_HISTONE_LABELS
#print(ks_results)
write.csv(ks_results, EPIGENETIC_MEAN_KS_STAT_SAMPLE_FILE)


#########################################################
# Select data for analysis and run Spearman correlation #
allrna_corr_obj <- rcorr(
  as.matrix(
    data_numeric_epigenetic_sample %>%
      #filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3",
      #                      "lncrna-exon1", "lncrna-exon2",
      #                      "short-ncrna")) %>%
      select(all_of(SELECTED_HISTONE_FEATURES))
  ),
  type = "spearman"
)

# Unpack correlation coefficients and p-values #
allrna_corr_matrix <- allrna_corr_obj[["r"]]
allrna_corr_matrix[is.na(allrna_corr_matrix)] <- NA
allrna_pvalues <- allrna_corr_obj[["P"]]

# Number of features
n <- ncol(allrna_corr_matrix)
# Number of unique tests
number_of_tests <- n * (n - 1) / 2
# Original alpha
original_alpha <- 0.05
# Corrected alpha
corrected_alpha <- original_alpha / number_of_tests

# Adjust p-values using Benjamini-Hochberg
#help("p.adjust")
allrna_p_values_adjusted <- p.adjust(allrna_pvalues, method="BH")

allrna_p_values_adjusted <- matrix(allrna_p_values_adjusted, nrow = n, byrow = TRUE)
diag= rep(0,n)
diag(allrna_p_values_adjusted) <- diag

# Use short names
rownames(allrna_corr_matrix) <- SELECTED_HISTONE_LABELS
colnames(allrna_corr_matrix) <- SELECTED_HISTONE_LABELS

# Use short names
rownames(allrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS
colnames(allrna_p_values_adjusted) <- SELECTED_HISTONE_LABELS


#########################
# Use heatmaply package #
# Calculate hierarchical clustering for rows and columns
hc_rows <- hclust(dist(abs(allrna_corr_matrix))) # Compute absolute value for clustering.
hc_cols <- hclust(dist(t(abs(allrna_corr_matrix)))) # Compute absolute value for clustering.

# Create dendrograms
row_dendro <- as.dendrogram(hc_rows)
col_dendro <- as.dendrogram(hc_cols)

# Get the column order from the clustering to align our plots later
col_order <- order.dendrogram(col_dendro)
row_order <- order.dendrogram(row_dendro)

# 3. Create the matrix for cell notes
cellnote_matrix <- matrix(NA_character_, nrow = nrow(allrna_corr_matrix), ncol = ncol(allrna_corr_matrix))
if (!is.null(dimnames(allrna_corr_matrix))) {
  dimnames(cellnote_matrix) <- dimnames(allrna_corr_matrix)
}

for (r_idx in 1:nrow(allrna_corr_matrix)) {
  for (c_idx in 1:ncol(allrna_corr_matrix)) {
    corr_val <- allrna_corr_matrix[r_idx, c_idx]
    p_val_adj <- allrna_p_values_adjusted[r_idx, c_idx]
    
    if (!is.na(corr_val)) { # This is the key condition
      text_val <- sprintf("%.1f", corr_val) 
      
      if (!is.na(p_val_adj) && p_val_adj >= corrected_alpha) {
        cellnote_matrix[r_idx, c_idx] <- " X"
      } else {
        cellnote_matrix[r_idx, c_idx] <- paste0("",text_val)
      }
    } else {
      # If corr_val is NA, cellnote_matrix entry remains NA_character_ (its initialized value)
      # cellnote_matrix[r_idx, c_idx] <- NA_character_ # Explicitly, though already NA
    }
  }
}

# Create the heatmap with dendrograms using heatmaply
heatmaply(allrna_corr_matrix, 
          Rowv = row_dendro, 
          Colv = col_dendro, 
          #col_side_colors = col_annotation_matrix,
          cellnote = cellnote_matrix, 
          cellnote_size = 20, 
          symm = TRUE,
          colors = colorRampPalette(c("navy", "white", "red"))(50), 
          limits = c(-1, 1), # MODIFIED: Explicitly set color scale limits
          grid_color = "grey80", 
          cellnote_textposition = "middle center", 
          grid_width = 0.001,    
          fontsize_row = 20,
          fontsize_col = 20,
          margins = c(60, 60, 40, 20) 
          )


#####################################
#####################################
# VERSION 4, Add ks stats at the bottom 
# 3. Create the TWO individual plot components

# Plot A: Main Heatmap (Center)
main_heatmap <- heatmaply(
  allrna_corr_matrix,
  Rowv = row_dendro, Colv = col_dendro,
  show_dendrogram = c(FALSE, TRUE), 
  hide_colorbar = TRUE,
  
  cellnote = cellnote_matrix,
  cellnote_size = 20, cellnote_textposition = "middle center",
  
  symm = TRUE,
  colors = colorRampPalette(c("navy", "white", "red"))(50),
  limits = c(-1, 1),
  grid_color = "grey80", grid_width = 0.001,
  
  fontsize_row = 20, fontsize_col = 20,
  margins = c(0, 0, 0, 0)
)

# Plot B: K-S D value annotation heatmap (Bottom)
ks_d_matrix <- data.frame(
  "mRNA KS" = ks_d_values_m$signed_D,
  "sncRNA KS" = ks_d_values_s$signed_D,
  "lncRNA KS" = ks_d_values_l$signed_D,
  "mean|KS|" = ks_results$meanAbsD,
  check.names = FALSE)

ks_d_matrix <- t(ks_d_matrix)
colnames(ks_d_matrix) <- SELECTED_HISTONE_LABELS

# Create cellnotes for the new multi-row matrix
ks_cellnote <- matrix(sprintf("%.2f", ks_d_matrix), nrow = nrow(ks_d_matrix), dimnames = dimnames(ks_d_matrix))

# We must reorder the columns to match the main heatmap's clustering
ks_d_matrix_ordered <- ks_d_matrix[, col_order, drop = FALSE]
ks_cellnote_ordered <- ks_cellnote[, col_order, drop = FALSE]

ks_heatmap <- heatmaply(
  ks_d_matrix_ordered,
  Rowv = FALSE, Colv = FALSE,
  cellnote = ks_cellnote_ordered,
  cellnote_color = "black", cellnote_size = 20, cellnote_textposition = "middle center",
  
  colors = colorRampPalette(c("navy", "white", "red"))(50),
  limits = c(-1, 1), # Use the same limits as the main heatmap for consistent color scale
  grid_color = "grey80", grid_width = 0.001,
  
  dendrogram = "none",
  # No longer using a single ylab, the row names of the matrix will be used
  
  fontsize_row = 20, fontsize_col = 20,
  margins = c(0, 0, 0, 0) # Margin on bottom for column labels
)


# 4. Combine the plots into a final layout
# Stack the two plots.
right_stack <- subplot(main_heatmap, ks_heatmap, 
                       nrows = 2, 
                       shareX = TRUE, 
                       margin = 0.005,
                       heights = c(0.9, 0.1), # Dendrogram, Main, Annotation
                       titleY = TRUE)


# Display the final plot interactively
right_stack

# Save static exports via webshot2 (requires Chrome; no Python/kaleido needed)
tmp_html <- tempfile(fileext = ".html")
saveWidget(right_stack, tmp_html, selfcontained = FALSE)

options(chromote.timeout = 60)

webshot2::webshot(tmp_html,
                  file  = file.path(VIOLIN_PLOTS_DIR, "png", "figureS3-epigenetic-correlation-matrix.png"),
                  vwidth = 3840, vheight = 2160, delay = 5)

# Use chromote directly: set viewport to match PNG, then print to exact-size single-page PDF
ch <- chromote::ChromoteSession$new()
ch$Emulation$setDeviceMetricsOverride(
  width = 3840, height = 2160, deviceScaleFactor = 1, mobile = FALSE
)
ch$Page$navigate(paste0("file:///", gsub("\\\\", "/", normalizePath(tmp_html))))
Sys.sleep(6)
pdf_data <- ch$Page$printToPDF(
  printBackground = TRUE,
  paperWidth      = 40,    # 3840px / 96dpi = 40 inches wide  (landscape, no rotation needed)
  paperHeight     = 22.5,  # 2160px / 96dpi = 22.5 inches tall
  marginTop       = 0, marginBottom = 0, marginLeft = 0, marginRight = 0
)
writeBin(jsonlite::base64_dec(pdf_data$data),
         file.path(VIOLIN_PLOTS_DIR, "pdf", "figureS3-epigenetic-correlation-matrix.pdf"))
ch$close()

unlink(tmp_html)
