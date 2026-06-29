# Violin + boxplot of SNP density per RNA type.
# Modernized to run against the current feature pipeline: load_gene_functionality_features()
# now returns a single combined feature_matrix (with a Dataset column) rather than globals,
# and the rna_type annotation lives at data/rnatype_feature/. Mirrors scripts/daf_vs_rnatype_plot.R.
#
# Extended with three pooled functional-sncRNA groups in addition to the per-type violins:
#   - "Pol II"        : functional sncRNA transcribed by RNA Pol II
#   - "Pol III"       : functional sncRNA transcribed by RNA Pol III
#   - "all sncRNA(+)" : all functional sncRNA pooled (positive counterpart to sncRNA(-))
# pol_class comes from data/daf/sncrna_positives.tsv (the only place it exists).
#
# Run from working_scripts/ :  Rscript scripts/snp_density_vs_rnatype_plot.R

library(dplyr)
library(ggplot2)

source("complementary/config.R")
source("complementary/load_gene_functionality_features.R")

feature_matrix <- load_gene_functionality_features()

# --- Functional sncRNA: attach rna_type (annotation) and pol_class (DAF positives) ----------
# The functional-sncRNA feature rows, the rna_type annotation CSV, and sncrna_positives.tsv are
# all 1000 rows in identical order (the project-wide positional alignment that
# daf_zscore_by_pol_class_plot.R also relies on). Guard it with stopifnot().
sncrna_func <- feature_matrix %>% dplyr::filter(Dataset == "short-ncrna")
sncrna_ann  <- read.csv("complementary/rnatype_feature/functional-short-ncrna-rnatype-hgncid-rnacid-feature.csv",
                        stringsAsFactors = FALSE)
sncrna_pos  <- read.delim(DAF_SNCRNA_POS_FILE, stringsAsFactors = FALSE)
stopifnot(nrow(sncrna_func) == nrow(sncrna_ann),
          nrow(sncrna_func) == nrow(sncrna_pos))
sncrna_func$rna_type  <- sncrna_ann$rna_type
sncrna_func$pol_class <- sncrna_pos$pol_class

# --- Assemble per-type rows (SNP_density + rna_type) ----------------------------------------
mrna <- feature_matrix %>%
  dplyr::filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>%
  dplyr::transmute(SNP_density, rna_type = "mRNA")
lncrna <- feature_matrix %>%
  dplyr::filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")) %>%
  dplyr::transmute(SNP_density, rna_type = "lncRNA")
sncrna_nc <- feature_matrix %>%
  dplyr::filter(Dataset == "short-ncrna-negative-control") %>%
  dplyr::transmute(SNP_density, rna_type = "sncRNA(-)")
sncrna_types <- sncrna_func %>% dplyr::transmute(SNP_density, rna_type)

# --- Three new pooled functional-sncRNA groups ----------------------------------------------
pol2_pool <- sncrna_func %>% dplyr::filter(pol_class == "Pol II")  %>%
  dplyr::transmute(SNP_density, rna_type = "sncRNA Pol II")
pol3_pool <- sncrna_func %>% dplyr::filter(pol_class == "Pol III") %>%
  dplyr::transmute(SNP_density, rna_type = "sncRNA Pol III")
allsnc_pool <- sncrna_func %>% dplyr::transmute(SNP_density, rna_type = "sncRNA(+)")

allDatasets <- dplyr::bind_rows(mrna, lncrna, sncrna_types, sncrna_nc,
                                pol2_pool, pol3_pool, allsnc_pool) %>%
  dplyr::filter(!is.na(SNP_density))

# Pretty labels for sncRNA subtypes (the annotation labels snaR genes "ncRNA").
allDatasets$rna_type <- dplyr::recode(allDatasets$rna_type,
  "pre_miRNA" = "pre-miRNA",
  "vault_RNA" = "vault RNA",
  "snaR_RNA"  = "snaR RNA",
  "ncRNA"     = "snaR RNA"
)

cat("Row counts per rna_type:\n")
print(table(allDatasets$rna_type, useNA = "ifany"))

# Order rna_type by median SNP density, ascending (existing behaviour).
medians <- allDatasets %>%
  dplyr::group_by(rna_type) %>%
  dplyr::summarise(median_density = median(SNP_density, na.rm = TRUE)) %>%
  dplyr::arrange(median_density)
allDatasets$rna_type <- factor(allDatasets$rna_type, levels = medians$rna_type)

# Project palette (CLAUDE.md §3); interpolated when categories exceed the base colours.
palette_full <- c("#F4A582FF", "#FCF2F1", "#FAECEA", "#e37b88FF", "#56bdfcFF",
                  "#F8E6E3", "#F7DFDC", "#F5D9D4", "#EBB0A6", "#E4988B")
n_levels <- nlevels(allDatasets$rna_type)
palette <- if (n_levels <= length(palette_full)) {
  palette_full[seq_len(n_levels)]
} else {
  grDevices::colorRampPalette(palette_full)(n_levels)
}

# Name the interpolated palette by category, then pin specific violins to the
# project-wide colours so the scheme matches the other violin plots. snoRNA takes
# the colour sncRNA(+) used before; sncRNA(+) becomes the canonical red and
# sncRNA(-) the canonical blue. All other categories keep their interpolated colours.
names(palette) <- levels(allDatasets$rna_type)
palette["snoRNA"]    <- "#F7E4E1"     # previous sncRNA(+) colour
palette["sncRNA(+)"] <- "#D6603FFF"
palette["sncRNA(-)"] <- "#56bdfcFF"

ggplot(allDatasets, aes(x = rna_type, y = SNP_density, fill = rna_type)) +
  geom_violin(scale = "width") +
  geom_boxplot(alpha = 0.0, outliers = FALSE,
               position = position_dodge(width = 0.9), width = 0.2, size = 1) +
  labs(title = "SNP density Distribution by RNA type", y = "SNP Density") +
  theme_minimal(base_size = 56) +
  theme(
    axis.text.x = element_text(size = 30, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 0),
    axis.title.y = element_text(size = 46),
    legend.position = "none",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 52, hjust = 0.5, margin = margin(b = 60)),
    plot.subtitle = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text = element_blank()
  ) +
  scale_fill_manual(values = palette) +
  coord_cartesian(ylim = c(0, 3)) # Zoom into a specific area without removing data points


ggsave("figure7-Violin-population-RNAType_Vs_SNPs.png", path = file.path(VIOLIN_PLOTS_DIR, "png"),
       scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
ggsave("figure7-Violin-population-RNAType_Vs_SNPs.pdf", path = file.path(VIOLIN_PLOTS_DIR, "pdf"),
       scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)

# --- KS tests: sncRNA(-) vs sncRNA Pol II / Pol III ------------------------------------
nc_snp   <- na.omit(sncrna_nc$SNP_density)
pol2_snp <- na.omit(pol2_pool$SNP_density)
pol3_snp <- na.omit(pol3_pool$SNP_density)

ks_output <- capture.output({
  cat(sprintf("Group sizes:  sncRNA(-) n=%d | Pol II n=%d | Pol III n=%d\n",
              length(nc_snp), length(pol2_snp), length(pol3_snp)))
  cat(sprintf("Medians:  sncRNA(-) %.4f | Pol II %.4f | Pol III %.4f\n\n",
              median(nc_snp), median(pol2_snp), median(pol3_snp)))
  cat("=== KS test: sncRNA(-) vs sncRNA Pol II ===\n")
  print(ks.test(nc_snp, pol2_snp))
  cat("\n=== KS test: sncRNA(-) vs sncRNA Pol III ===\n")
  print(ks.test(nc_snp, pol3_snp))
})

writeLines(ks_output, file.path(KS_STATS_DIR, "ks_snp_density_sncrna_neg_vs_pol_class.txt"))
cat(paste(ks_output, collapse = "\n"), "\n")
