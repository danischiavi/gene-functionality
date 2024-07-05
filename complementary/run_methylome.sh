#!/bin/bash

set -uex


# LONG NCRNA
## Process methylome marks for exon1
sh methylome_processing.sh ../data/functional-lncrna-exon1-dataset.csv lncrna-exon1-methylome-feature
sh methylome_processing.sh ../data/lncrna-exon1-negative-control-dataset-curated.csv lncrna-exon1-NC-methylome-feature
## Process methylome marks for exon2
sh methylome_processing.sh ../data/functional-lncrna-exon2-dataset.csv lncrna-exon2-methylome-feature
sh methylome_processing.sh ../data/lncrna-exon2-negative-control-dataset-curated.csv lncrna-exon2-NC-methylome-feature
## Join datasets into one for R script
sh join_datasets.sh ../data/methylome_feature/lncrna-exon1-methylome-feature.csv ../data/methylome_feature/lncrna-exon2-methylome-feature.csv ../data/methylome_feature/lncrna_positives_matrix
sh join_datasets.sh ../data/methylome_feature/lncrna-exon1-NC-methylome-feature.csv ../data/methylome_feature/lncrna-exon2-NC-methylome-feature.csv ../data/methylome_feature/lncrna_negatives_matrix
sh join_datasets.sh ../data/methylome_feature/lncrna_positives_matrix.csv ../data/methylome_feature/lncrna_negatives_matrix.csv ../data/methylome_feature/lncrna_matrix

#SHORT NCRNA
sh methylome_processing.sh \
    ../data/functional-short-ncrna-dataset.csv \
    short-ncrna-methylome-feature
sh methylome_processing.sh \
    ../data/short-ncrna-negative-control-dataset-curated.csv \
    short-ncrna-NC-methylome-feature
sh join_datasets.sh \
    ../data/methylome_feature/short-ncrna-methylome-feature.csv \
    ../data/methylome_feature/short-ncrna-NC-methylome-feature.csv \
    ../data/methylome_feature/short_ncrna_matrix


# PROTEIN CODING
sh methylome_processing.sh ../data/functional-protein-exon2-dataset.csv protein-exon2-methylome-feature
sh methylome_processing.sh ../data/protein-exon2-negative-control-dataset-curated.csv protein-exon2-NC-methylome-feature

sh methylome_processing.sh ../data/functional-protein-exon3-dataset.csv protein-exon3-methylome-feature
sh methylome_processing.sh ../data/protein-exon3-negative-control-dataset-curated.csv protein-exon3-NC-methylome-feature

sh join_datasets.sh ../data/methylome_feature/protein-exon2-methylome-feature.csv ../data/methylome_feature/protein-exon3-methylome-feature.csv ../data/methylome_feature/protein_positives_matrix
sh join_datasets.sh ../data/methylome_feature/protein-exon2-NC-methylome-feature.csv ../data/methylome_feature/protein-exon3-NC-methylome-feature.csv ../data/methylome_feature/protein_negatives_matrix
sh join_datasets.sh ../data/methylome_feature/protein_positives_matrix.csv ../data/methylome_feature/protein_negatives_matrix.csv ../data/methylome_feature/protein_matrix

