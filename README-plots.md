# Paper Figures — Data Source Mapping

---

## Table A — Main figures (`figure1`–`figure8`)

| Figure | Source | Columns used | Notes |
|---|---|---|---|
| `figure1-KS-across-features-and-dataset` | `results/*-features_20260608.csv` — all 10 functional + negative-control files | KS distance computed across all 31 feature columns, all 5 RNA-type/exon datasets | Generated with `scripts/ks-plots.py` |
| `figure2-Violin-transcriptome` | `results/z-scores/*-zscores.csv` — all 10 files | `RPKM_tissue`, `RPKM_primary.cell` | Generated with `complementary/generate_gene_functionality_zscore_violin_plots.R`. |
| `figure3-Violin-conservation` | `results/z-scores/*-zscores.csv` — all 10 files | `phyloP_max_241w`, `phyloP_max_100w`, `GERP_91_mammals_max`, `GERP_63_amniotes_max` | Generated with `complementary/generate_gene_functionality_zscore_violin_plots.R`. |
| `figure4-Violin-epigenetic` | `results/z-scores/*-zscores.csv` — all 10 files | `H3K9ac_MaxScaledSignal`, `H3K79me2_MaxScaledSignal`, `H3K79me1_MaxScaledSignal`, `chrm_acc_MaxScaledSignal`, `methylome` | Generated with `complementary/generate_gene_functionality_zscore_violin_plots.R`. |
| `figure5-Violin-RNAspecific` | `results/z-scores/*-zscores.csv` — all 10 files | `MFE`, `MaxEnt_splicing`, `fickett`, `coding_potential`, `Interaction_ave`, `Max_covariance` | Generated with `complementary/generate_gene_functionality_zscore_violin_plots.R`. |
| `figure6-Violin-repeats` | `results/z-scores/*-zscores.csv` — all 10 files | `copy_number`, `repeat_distance` | Generated with `complementary/generate_gene_functionality_zscore_violin_plots.R`. |
| `figure7-Violin-population` | `results/z-scores/*-zscores.csv` — all 10 files + `complementary/daf/sncrna_positives.tsv` + `complementary/rnatype_feature/functional-short-ncrna-rnatype-hgncid-rnacid-feature.csv` | `SNP_density`, `MAF_avg` | Generated with `complementary/generate_gene_functionality_zscore_violin_plots.R` + `complementary/snp_density_vs_rnatype_plot.R`. |
| `figure8-Violin-intrinsic` | `results/z-scores/*-zscores.csv` — all 10 files | `GC_percentage`, `GA`, `GT`, `CpG`, `TA`, `T`, `C`, `G`, `lowComplexity_density` | Generated with `complementary/generate_gene_functionality_zscore_violin_plots.R`. |

---

## Table B — Supplementary figures (`figureS1`, `figureS3`–`S13`)

| Figure | Source | Columns used | Notes |
|---|---|---|---|
| `figureS1-Upset-plot` | `data/raw/gencodev44_chr22.bed`,`data/raw/gencodev45_chr22.bed`,`data/raw/hgnc_chr22.bed`,`data/raw/ncbi_chr22.bed`,`data/raw/rnacentral_chr22.bed`,`data/raw/uniprot_chr22.bed` | | Generated with `complementary/UpsetPlot.R`. |
| `figureS3-epigenetic-correlation-matrix` | `complementary/histone_marks/*.csv` | `H*_MaxScaledSignal`, `chrm_acc_MaxScaledSignal`, `methylome` | Generated with `complementary/generate_epigenetic_heatmap.R`. |
| `figureS4-mRNA-correlation-matrix` | `results/functional-protein-exon2/exon3-...-features_*.csv` (+ negative controls) | | Feature-feature Spearman correlation for protein-coding, both exons combined. Generated with `complementary/generate_gene_functionality_heatmap.R`. |
| `figureS5-sncRNA-correlation-matrix` | `results/functional-short-ncrna-dataset-features_*.csv` + `results/short-ncrna-negative-control-dataset-features_*.csv` |  | Feature-feature correlation for short ncRNA. Generated with `complementary/generate_gene_functionality_heatmap.R`. |
| `figureS6-lncRNA-correlation-matrix` | `results/functional-lncrna-exon1/exon2-...-features_*.csv` (+ negative controls) |  | Feature-feature correlation for lncRNA, both exons combined. Generated with `complementarygenerate_gene_functionality_heatmap.R`. |
| `figureS7-data-effect-sizes` | `results/*-features_*.csv` — all 10 files |  | Effect size of functional vs. negative-control per feature, all RNA types. Generated with `scripts/ks-plots.py` 
| `figureS8-sncRNA-primary-cell-RPKM-z-score` | `results/z-scores/functional-short-ncrna-dataset-zscores.csv` + `results/z-scores/short-ncrna-negative-control-dataset-zscores.csv` | Column `RPKM_primary.cell`. | Generated with `complementary/generate_gene_functionality_zscore_distribution_plots.R`. |
| `figureS9-mRNA-PhyloP-mammals-z-score` | `results/z-scores/functional-protein-exon2/exon3-dataset-zscores.csv` (+ negative controls) | Column `phyloP_max_241w`. | Generated with `complementary/generate_gene_functionality_zscore_distribution_plots.R`. |
| `figureS10-mRNA-methylome-z-score` | `results/z-scores/functional-protein-exon2/exon3-dataset-zscores.csv` (+ negative controls) | Column `methylome`. | Generated with `complementary/generate_gene_functionality_zscore_distribution_plots.R`. |
| `figureS11-mRNA-copies-z-score` | `results/z-scores/functional-protein-exon2/exon3-dataset-zscores.csv` (+ negative controls) | Column `copy_number`. | Generated with `complementary/generate_gene_functionality_zscore_distribution_plots.R`. |
| `figureS12-lncRNA-copies-z-score` | `results/z-scores/functional-lncrna-exon1/exon2-dataset-zscores.csv` (+ negative controls) | Column `copy_number`. | Generated with `complementary/generate_gene_functionality_zscore_distribution_plots.R`. |
| `figureS13-Heatscatter-effect-distance` | Effect sizes from `results/*-features_*.csv` | `DistanceGene`, `GC%`, `SNP_density` | Generated with `complementary/distance_effect_analysis.R`. |

---

## Table C — `results/Paper figures/Distribution_plots/` (pdfs + pngs, one row per feature × lncRNA/mRNA/sncRNA)

Source for every row: `results/z-scores/functional-{type}-...-dataset_zscores.csv` + `results/z-scores/-{type}--negative-control-dataset-zscores.csv` (For lncRNA use exon1+exon2 combined; for mRNA use protein-exon2+exon3 combined). Generated with `complementary/generate_gene_functionality_zscore_distribution_plots.R`.

| Filename feature | Column | Filename feature | Column |
|---|---|---|---|
| `C` | `C` | `MAF` | `MAF_avg` |
| `CpG` | `CpG` | `MFE` | `MFE` |
| `G` | `G` | `MaxEnt_splicing` | `MaxEnt_splicing` |
| `GA` | `GA` | `RNAcode` | `coding_potential` |
| `GC%` | `GC_percentage` | `RPKM_primary_cell` | `RPKM_primary.cell` |
| `GERP_mammals` | `GERP_91_mammals_max` | `RPKM_tissue` | `RPKM_tissue` |
| `GERP_vertebrates` | `GERP_63_amniotes_max` | `Random` | `Random` |
| `GT` | `GT` | `SNPs` | `SNP_density` |
| `H3K79me1` | `H3K79me1_MaxScaledSignal` | `T` | `T` |
| `H3K79me2` | `H3K79me2_MaxScaledSignal` | `TA` | `TA` |
| `H3K9ac` | `H3K9ac_MaxScaledSignal` | `chrm_acc` | `chrm_acc_MaxScaledSignal` |
| `Interactions` | `Interaction_ave` | `copies` | `copy_number` |
| `covariance` | `Max_covariance` | `fickett` | `fickett` |
| `lowComplexity` | `lowComplexity_density` | `methylome` | `methylome` |
| `phyloP_mammals` | `phyloP_max_241w` | `phyloP_vertebrates` | `phyloP_max_100w` |
| `repeat_free` | `repeat_distance` | | |

---