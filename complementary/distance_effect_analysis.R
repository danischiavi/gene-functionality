# Effect on Distance Analysis
# Try to identify effects on distance for features of functional genes
# Load libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(patchwork)
library(LSD)

source("complementary/config.R")
source("complementary/load_gene_functionality_zscores.R")

dir.create(LATEST_LOG_SCATTER_PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
zscores_all <- load_gene_functionality_zscores()

# Check the column names of the loaded data
colnames(zscores_all)

# List features of interest (in this case, the 28 most informative)
eod_select_features <- c("GC_percentage", "lowComplexity_density",
                         "CpG", "GA", "GT", "TA",
                         "T","C","G",
                         "phyloP_max_241w", "phyloP_max_100w", 
                         "GERP_91_mammals_max", "GERP_63_amniotes_max", 
                         "RPKM_tissue", "RPKM_primary.cell", 
                         "H3K9ac_MaxScaledSignal", "H3K79me2_MaxScaledSignal", "H3K79me1_MaxScaledSignal", 
                         "chrm_acc_MaxScaledSignal", "methylome", 
                         "repeat_distance", "copy_number", "coding_potential", 
                         "fickett", "Max_covariance", "MFE", "MaxEnt_splicing",
                         "Interaction_ave", "SNP_density", "MAF_avg"
                         )

eod_select_features_labels <- c("GC%", "Complexity",
                                "CpG", "GA", "GT", "TA",
                                "T%","C%","G%",
                                "PhyloP-mammals", "PhyloP-vertebrates",
                                "GERP-mammals", "GERP-vertebrates",
                                "Tissue RPKM", "Primary cell RPKM", 
                                "H3K9ac", "H3K79me2", "H3K79me1",
                                "Chromatin", "Methylome",
                                "Repeat.free", "Copies", "RNAcode",
                                "Fickett", "Covariance", "MFE", "MaxEnt splicing",
                                "Interactions", "SNPs", "MAF"
                                )


# Select just desired features
unique(zscores_all$Dataset)
funcProtExon2DataSelect <- zscores_all %>% filter(Dataset == "protein-coding-exon2") %>% dplyr::select(all_of(eod_select_features))
funcProtExon3DataSelect <- zscores_all %>% filter(Dataset == "protein-coding-exon3") %>% dplyr::select(all_of(eod_select_features))
funcLncrnaExon1DataSelect <- zscores_all %>% filter(Dataset == "lncrna-exon1") %>% dplyr::select(all_of(eod_select_features))
funcLncrnaExon2DataSelect <- zscores_all %>% filter(Dataset == "lncrna-exon2") %>% dplyr::select(all_of(eod_select_features))
funcSncrnaDatasetSelect <- zscores_all %>% filter(Dataset == "short-ncrna") %>% dplyr::select(all_of(eod_select_features))

protExon2NCDataSelect <- zscores_all %>% filter(Dataset == "protein-exon2-negative-control") %>% dplyr::select(all_of(eod_select_features))
protExon3NCDataSelect <- zscores_all %>% filter(Dataset == "protein-exon3-negative-control") %>% dplyr::select(all_of(eod_select_features))
lncrnaExon1NCDataSelect <- zscores_all %>% filter(Dataset == "lncrna-exon1-negative-control") %>% dplyr::select(all_of(eod_select_features))
lncrnaExon2NCDataSelect <- zscores_all %>% filter(Dataset == "lncrna-exon2-negative-control") %>% dplyr::select(all_of(eod_select_features))
sncrnaNCDataSelect <- zscores_all %>% filter(Dataset == "short-ncrna-negative-control") %>% dplyr::select(all_of(eod_select_features))


# Retrieve DistanceGene feature and add to negative control z-score dataframes
protExon2NCDataSelect$DistanceGene <- as.data.frame(read.csv(NC_PROT_EXON2_FEATURES_FILE, header = TRUE))$DistanceGene
protExon3NCDataSelect$DistanceGene <- as.data.frame(read.csv(NC_PROT_EXON3_FEATURES_FILE, header = TRUE))$DistanceGene
lncrnaExon1NCDataSelect$DistanceGene <- as.data.frame(read.csv(NC_LNCRNA_EXON1_FEATURES_FILE, header = TRUE))$DistanceGene
lncrnaExon2NCDataSelect$DistanceGene <- as.data.frame(read.csv(NC_LNCRNA_EXON2_FEATURES_FILE, header = TRUE))$DistanceGene
sncrnaNCDataSelect$DistanceGene <- as.data.frame(read.csv(NC_SNCRNA_FEATURES_FILE, header = TRUE))$DistanceGene


# Create scatter plots with LOESS curve #
# PROTEIN: #
prot_pos_zscores <- zscores_all %>%
  filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>%
  mutate(distance=1)

prot_neg_zscores <- zscores_all %>%
  filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>%
  mutate(distance = rbind(
    protExon2NCDataSelect %>% dplyr::select(DistanceGene),
    protExon3NCDataSelect %>% dplyr::select(DistanceGene)
  )$DistanceGene)


prot_pos_zscores$group <- "mrna(+)"
prot_neg_zscores$group <- "mrna(-)"
prot_zscores_all <- rbind(prot_pos_zscores,prot_neg_zscores)
prot_zscores_all$group <- factor(prot_zscores_all$group, levels = c("mrna(+)","mrna(-)"))
table(prot_zscores_all$group)
# GC%
prot_gc_plot <- ggplot(prot_neg_zscores, aes(x = distance, y = GC_percentage, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "mRNA") +
  xlab("Distance to Gene (Mb)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    #legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  scale_x_continuous(trans='log10', breaks = c(1000, 10000, 100000, 1000000, 5000000), 
                     labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"),
                     limits = c(1000, 5000000)) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c("dodgerblue", "orange"))
prot_gc_plot
#ggsave(paste0("distance_effect_protein_loess_wiggle_veplus",'_',"GCcontent",'.png'), path = LATEST_LOG_SCATTER_PLOTS_DIR, prot_gc_plot, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# SNP density
prot_snp_plot <- ggplot(prot_zscores_all, aes(x = log10(distance), y = SNP_density, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "mRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,4)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
prot_snp_plot
#ggsave(paste0("distance_effect_protein_loess_wiggle_veplus",'_',"SNPdensity",'.png'), path = LATEST_LOG_SCATTER_PLOTS_DIR, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)



# LNCRNA: #
lncrna_pos_zscores <- zscores_all %>%
  filter(Dataset == "lncrna-exon1" | Dataset == "lncrna-exon2") %>%
  mutate(distance=1)

lncrna_neg_zscores <- zscores_all %>%
  filter(Dataset == "lncrna-exon1-negative-control" | Dataset == "lncrna-exon2-negative-control") %>%
  mutate(distance = rbind(
    lncrnaExon1NCDataSelect %>% dplyr::select(DistanceGene),
    lncrnaExon2NCDataSelect %>% dplyr::select(DistanceGene)
  )$DistanceGene)

lncrna_pos_zscores$group <- "lncrna(+)"
lncrna_neg_zscores$group <- "lncrna(-)"
lncrna_zscores_all <- rbind(lncrna_pos_zscores,lncrna_neg_zscores)
lncrna_zscores_all$group <- factor(lncrna_zscores_all$group, levels = c("lncrna(+)","lncrna(-)"))

# GC%
lncrna_gc_plot <- ggplot(lncrna_zscores_all, aes(x = log10(distance), y = GC_percentage, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "lncRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
lncrna_gc_plot
#ggsave(paste0("distance_effect_lncrna_loess_wiggle_veplus",'_',"GCcontent",'.png'), path = LATEST_LOG_SCATTER_PLOTS_DIR, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


# SNP density
lncrna_snp_plot <- ggplot(lncrna_zscores_all, aes(x = log10(distance), y = SNP_density, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "lncRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,4)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
lncrna_snp_plot
#ggsave(paste0("distance_effect_lncrna_loess_wiggle_veplus",'_',"SNPdensity",'.png'), path = LATEST_LOG_SCATTER_PLOTS_DIR, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


# SNCRNA: #
sncrna_pos_zscores <- zscores_all %>%
  filter(Dataset == "short-ncrna") %>%
  mutate(distance=1)

sncrna_neg_zscores <- zscores_all %>%
  filter(Dataset == "short-ncrna-negative-control") %>%
  mutate(distance = (sncrnaNCDataSelect %>% 
           dplyr::select(DistanceGene))$DistanceGene
  )


sncrna_pos_zscores$group <- "sncrna(+)"
sncrna_neg_zscores$group <- "sncrna(-)"
sncrna_zscores_all <- rbind(sncrna_pos_zscores,sncrna_neg_zscores)
sncrna_zscores_all$group <- factor(sncrna_zscores_all$group, levels = c("sncrna(+)","sncrna(-)"))

# GC%
sncrna_gc_plot <- ggplot(sncrna_zscores_all, aes(x = log10(distance), y = GC_percentage, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "sncRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
#ggsave(paste0("distance_effect_sncrna_loess_wiggle_veplus",'_',"GCcontent",'.png'), path = LATEST_LOG_SCATTER_PLOTS_DIR, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


# SNP density
sncrna_snp_plot <- ggplot(sncrna_zscores_all, aes(x = log10(distance), y = SNP_density, color = group)) +
  geom_point(show.legend = TRUE) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1, span = 0.5) +
  facet_wrap(~ "", scales = "free") +
  labs(subtitle = "sncRNA") +
  xlab("log(Distance to Gene)") +
  ylab("Robust Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 30, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 36),
    axis.title.y = element_text(size = 36),
    legend.position = "none",
    plot.subtitle = element_text(size = 42, hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-2,4)) +
  scale_color_manual(values = c("orange", "dodgerblue"))
#ggsave(paste0("distance_effect_sncrna_loess_wiggle_veplus",'_',"SNPdensity",'.png'), path = LATEST_LOG_SCATTER_PLOTS_DIR, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)


##############################
### USE heatscatter R PACKAGE.

# Wrap all drawing calls in a function so they can be replayed into both
# a png() device and a pdf() device without duplicating code.
# R's lexical scoping gives the function access to prot_neg_zscores etc.
draw_heatscatter_plots <- function() {
  layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))

  # GC%
  par(mar = c(7, 5, 7, 2))
  heatscatter(prot_neg_zscores$distance,
              prot_neg_zscores$GC_percentage,
              daltonize = TRUE, cor = FALSE, log = "x",
              xlim = c(7500, 5100000), ylim = c(-3, 3),
              xlab = "Distance to Gene", ylab = "Robust Z-score",
              main = "", xaxt = 'n', cex.axis = 2, cex.lab = 2.5)
  axis(1, at = c(1000, 10000, 100000, 1000000, 5000000),
       labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"), cex.axis = 2)
  title(main = "mRNA", cex.main = 3)
  valid_data <- !is.na(prot_neg_zscores$distance) & !is.na(prot_neg_zscores$GC_percentage)
  sorted_idx <- order(prot_neg_zscores$distance[valid_data])
  sorted_distance <- prot_neg_zscores$distance[valid_data][sorted_idx]
  sorted_gc <- prot_neg_zscores$GC_percentage[valid_data][sorted_idx]
  loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
  pred <- predict(loess_fit, se = TRUE)
  lines(sorted_distance, pred$fit, col = "red", lwd = 3)
  lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
  lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

  par(mar = c(7, 5, 7, 2))
  heatscatter(sncrna_neg_zscores$distance,
              sncrna_neg_zscores$GC_percentage,
              daltonize = TRUE, cor = FALSE, log = "x",
              xlim = c(7500, 5100000), ylim = c(-3, 3),
              xlab = "Distance to Gene", ylab = "Robust Z-score",
              main = "", xaxt = 'n', cex.axis = 2, cex.lab = 2.5)
  axis(1, at = c(1000, 10000, 100000, 1000000, 5000000),
       labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"), cex.axis = 2)
  title(main = "sncRNA", cex.main = 3)
  valid_data <- !is.na(sncrna_neg_zscores$distance) & !is.na(sncrna_neg_zscores$GC_percentage)
  sorted_idx <- order(sncrna_neg_zscores$distance[valid_data])
  sorted_distance <- sncrna_neg_zscores$distance[valid_data][sorted_idx]
  sorted_gc <- sncrna_neg_zscores$GC_percentage[valid_data][sorted_idx]
  loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
  pred <- predict(loess_fit, se = TRUE)
  lines(sorted_distance, pred$fit, col = "red", lwd = 3)
  lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
  lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

  heatscatter(lncrna_neg_zscores$distance,
              lncrna_neg_zscores$GC_percentage,
              daltonize = TRUE, cor = FALSE, log = "x",
              xlim = c(7500, 5100000), ylim = c(-3, 3),
              xlab = "Distance to Gene", ylab = "Robust Z-score",
              main = "", xaxt = 'n', cex.axis = 2, cex.lab = 2.5)
  axis(1, at = c(1000, 10000, 100000, 1000000, 5000000),
       labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"), cex.axis = 2)
  title(main = "lncRNA", cex.main = 3)
  valid_data <- !is.na(lncrna_neg_zscores$distance) & !is.na(lncrna_neg_zscores$GC_percentage)
  sorted_idx <- order(lncrna_neg_zscores$distance[valid_data])
  sorted_distance <- lncrna_neg_zscores$distance[valid_data][sorted_idx]
  sorted_gc <- lncrna_neg_zscores$GC_percentage[valid_data][sorted_idx]
  loess_fit <- loess(sorted_gc ~ sorted_distance, span = 0.1)
  pred <- predict(loess_fit, se = TRUE)
  lines(sorted_distance, pred$fit, col = "red", lwd = 3)
  lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
  lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

  # SNP density
  heatscatter(prot_neg_zscores$distance,
              prot_neg_zscores$SNP_density,
              daltonize = TRUE, cor = FALSE, log = "x",
              xlim = c(7500, 5100000), ylim = c(-3, 3),
              xlab = "Distance to Gene", ylab = "Robust Z-score",
              main = "", xaxt = 'n', cex.axis = 2, cex.lab = 2.5)
  axis(1, at = c(1000, 10000, 100000, 1000000, 5000000),
       labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"), cex.axis = 2)
  title(main = "mRNA", cex.main = 3)
  valid_data <- !is.na(prot_neg_zscores$distance) &
    !is.na(prot_neg_zscores$SNP_density) &
    prot_neg_zscores$SNP_density >= -3 &
    prot_neg_zscores$SNP_density <= 3
  sorted_idx <- order(prot_neg_zscores$distance[valid_data])
  sorted_distance <- prot_neg_zscores$distance[valid_data][sorted_idx]
  sorted_snp <- prot_neg_zscores$SNP_density[valid_data][sorted_idx]
  loess_fit <- loess(sorted_snp ~ sorted_distance, span = 0.1)
  pred <- predict(loess_fit, se = TRUE)
  lines(sorted_distance, pred$fit, col = "red", lwd = 3)
  lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
  lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

  heatscatter(sncrna_neg_zscores$distance,
              sncrna_neg_zscores$SNP_density,
              daltonize = TRUE, cor = FALSE, log = "x",
              xlim = c(7500, 5100000), ylim = c(-3, 3),
              xlab = "Distance to Gene", ylab = "Robust Z-score",
              main = "", xaxt = 'n', cex.axis = 2, cex.lab = 2.5)
  axis(1, at = c(1000, 10000, 100000, 1000000, 5000000),
       labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"), cex.axis = 2)
  title(main = "sncRNA", cex.main = 3)
  valid_data <- !is.na(sncrna_neg_zscores$distance) &
    !is.na(sncrna_neg_zscores$SNP_density) &
    sncrna_neg_zscores$SNP_density >= -3 &
    sncrna_neg_zscores$SNP_density <= 3
  sorted_idx <- order(sncrna_neg_zscores$distance[valid_data])
  sorted_distance <- sncrna_neg_zscores$distance[valid_data][sorted_idx]
  sorted_snp <- sncrna_neg_zscores$SNP_density[valid_data][sorted_idx]
  loess_fit <- loess(sorted_snp ~ sorted_distance, span = 0.1)
  pred <- predict(loess_fit, se = TRUE)
  lines(sorted_distance, pred$fit, col = "red", lwd = 3)
  lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
  lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)

  heatscatter(lncrna_neg_zscores$distance,
              lncrna_neg_zscores$SNP_density,
              daltonize = TRUE, cor = FALSE, log = "x",
              xlim = c(7500, 5100000), ylim = c(-3, 3),
              xlab = "Distance to Gene", ylab = "Robust Z-score",
              main = "", xaxt = 'n', cex.axis = 2, cex.lab = 2.5)
  axis(1, at = c(1000, 10000, 100000, 1000000, 5000000),
       labels = c("1kb", "10kb", "100kb", "1Mb", "5Mb"), cex.axis = 2)
  title(main = "lncRNA", cex.main = 3)
  valid_data <- !is.na(lncrna_neg_zscores$distance) &
    !is.na(lncrna_neg_zscores$SNP_density) &
    lncrna_neg_zscores$SNP_density >= -3 &
    lncrna_neg_zscores$SNP_density <= 3
  sorted_idx <- order(lncrna_neg_zscores$distance[valid_data])
  sorted_distance <- lncrna_neg_zscores$distance[valid_data][sorted_idx]
  sorted_snp <- lncrna_neg_zscores$SNP_density[valid_data][sorted_idx]
  loess_fit <- loess(sorted_snp ~ sorted_distance, span = 0.1)
  pred <- predict(loess_fit, se = TRUE)
  lines(sorted_distance, pred$fit, col = "red", lwd = 3)
  lines(sorted_distance, pred$fit + 1.96 * pred$se.fit, col = "grey", lty = 3)
  lines(sorted_distance, pred$fit - 1.96 * pred$se.fit, col = "grey", lty = 3)
}

# Close ALL open graphics devices (ggplot prints above may leave devices open)
graphics.off()

png(file.path(VIOLIN_PLOTS_DIR, "png", "figureS13-Heatscatter-effect-distance.png"),
    width=4920, height=3240, res=300)
draw_heatscatter_plots()
dev.off()

pdf(file.path(VIOLIN_PLOTS_DIR, "pdf", "figureS13-Heatscatter-effect-distance.pdf"),
    width=4920/300, height=3240/300)
draw_heatscatter_plots()
dev.off()

#################################
# Spearman Correlation Analysis #
# Function to compute a Spearman correlation coefficient for the specified features
compute_correlation <- function(zscores, selected_features) {
  corr_matrix_list <- list()
  for(feature in selected_features) {
    print(feature)
    
    corr_obj <- rcorr(zscores$distance,
                      zscores[[feature]], 
                      type = "spearman")
    corr_matrix_list[[feature]] <- list(rho = corr_obj[["r"]][2,1], p.value = corr_obj[["P"]][2,1])
  }
  return(corr_matrix_list)
}

# Function to Correct p-values
adjust_pvalues <- function(corr_matrix, selected_features) {
  pvalues <- list()
  for (feature in selected_features) {
    pvalues[[feature]] <- corr_matrix[[feature]]$p.value
  }
  pvalues_adj <- p.adjust(pvalues, method="BH", n = length(pvalues) * 3) # n equal length by 3 because we run a correlation test once per gene type (mRNA, sncRNA, and lncRNA)
  return(pvalues_adj)
}

# Function to Extract Spearman's rhos
extract_rhos <- function(corr_matrix, selected_features) {
  rhos <- list()
  for (feature in selected_features) {
    rhos[[feature]] <- corr_matrix[[feature]]$rho
  }
  return(rhos)
}

# Compute correlation coefficient to reveal distance effect:
#Protein coding
prot_corr_matrix <- compute_correlation(prot_neg_zscores, eod_select_features)
prot_rhos <- as.data.frame(extract_rhos(prot_corr_matrix, eod_select_features), check.names = FALSE)
rownames(prot_rhos) <- "rho"
prot_pvalues_adj <- as.data.frame(t(adjust_pvalues(prot_corr_matrix, eod_select_features)))
rownames(prot_pvalues_adj) <- "p.value"
prot_corr_matrix <- rbind(prot_rhos, prot_pvalues_adj)

#sncRNA
sncrna_corr_matrix <- compute_correlation(sncrna_neg_zscores, eod_select_features)
sncrna_rhos <- as.data.frame(extract_rhos(sncrna_corr_matrix, eod_select_features), check.names = FALSE)
rownames(sncrna_rhos) <- "rho"
sncrna_pvalues_adj <- as.data.frame(t(adjust_pvalues(sncrna_corr_matrix, eod_select_features)))
rownames(sncrna_pvalues_adj) <- "p.value"
sncrna_corr_matrix <- rbind(sncrna_rhos, sncrna_pvalues_adj)

#lncRNA
lncrna_corr_matrix <- compute_correlation(lncrna_neg_zscores, eod_select_features)
lncrna_rhos <- as.data.frame(extract_rhos(lncrna_corr_matrix, eod_select_features), check.names = FALSE)
rownames(lncrna_rhos) <- "rho"
lncrna_pvalues_adj <- as.data.frame(t(adjust_pvalues(lncrna_corr_matrix, eod_select_features)))
rownames(lncrna_pvalues_adj) <- "p.value"
lncrna_corr_matrix <- rbind(lncrna_rhos, lncrna_pvalues_adj)

prot_corr_matrix$measurement <- row.names(prot_corr_matrix)
sncrna_corr_matrix$measurement <- row.names(sncrna_corr_matrix)
lncrna_corr_matrix$measurement <- row.names(lncrna_corr_matrix)

prot_corr_matrix <- prot_corr_matrix |>
  pivot_longer(cols = "GC_percentage":"MAF_avg",
               names_to = "feature")

prot_corr_matrix <- prot_corr_matrix |>
  pivot_wider(names_from = measurement,
              values_from = value)

sncrna_corr_matrix <- sncrna_corr_matrix |>
  pivot_longer(cols = "GC_percentage":"MAF_avg",
               names_to = "feature")

sncrna_corr_matrix <- sncrna_corr_matrix |>
  pivot_wider(names_from = measurement,
              values_from = value)

lncrna_corr_matrix <- lncrna_corr_matrix |>
  pivot_longer(cols = "GC_percentage":"MAF_avg",
               names_to = "feature")

lncrna_corr_matrix <- lncrna_corr_matrix |>
  pivot_wider(names_from = measurement,
              values_from = value)

write.csv(prot_corr_matrix, file = "complementary/distance_effect_corr_matrix_protein.csv")
write.csv(sncrna_corr_matrix, file = "complementary/distance_effect_corr_matrix_sncrna.csv")
write.csv(lncrna_corr_matrix, file = "complementary/distance_effect_corr_matrix_lncrna.csv")
