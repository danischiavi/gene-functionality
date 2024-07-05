#!/bin/bash

set -uex

# Check for correct usage
if [ $# -ne 1 ]; then
    echo "Usage: $0 histone_name"
    exit 1
fi

# Store input file
HISTONE_NAME=$1


# LONG NCRNA
## Process histone marks for exon1
sh histone_processing.sh ../data/functional-lncrna-exon1-dataset.csv "$HISTONE_NAME" lncrna-exon1-histone-feature
sh histone_processing.sh ../data/lncrna-exon1-negative-control-dataset-curated.csv "$HISTONE_NAME" lncrna-exon1-NC-histone-feature
## Process histone marks for exon2
sh histone_processing.sh ../data/functional-lncrna-exon2-dataset.csv "$HISTONE_NAME" lncrna-exon2-histone-feature
sh histone_processing.sh ../data/lncrna-exon2-negative-control-dataset-curated.csv "$HISTONE_NAME" lncrna-exon2-NC-histone-feature
## Join datasets into one for R script
sh join_datasets.sh ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_lncrna-exon1-histone-feature.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_lncrna-exon2-histone-feature.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_lncrna_positives_matrix
sh join_datasets.sh ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_lncrna-exon1-NC-histone-feature.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_lncrna-exon2-NC-histone-feature.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_lncrna_negatives_matrix
sh join_datasets.sh ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_lncrna_positives_matrix.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_lncrna_negatives_matrix.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_lncrna_matrix

#SHORT NCRNA
sh histone_processing.sh \
    ../data/functional-short-ncrna-dataset.csv \
    "$HISTONE_NAME" \
    short-ncrna-histone-feature
sh histone_processing.sh \
    ../data/short-ncrna-negative-control-dataset-curated.csv \
    "$HISTONE_NAME" \
    short-ncrna-NC-histone-feature
sh join_datasets.sh \
    ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_short-ncrna-histone-feature.csv \
    ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_short-ncrna-NC-histone-feature.csv \
    ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_short_ncrna_matrix


# PROTEIN CODING
sh histone_processing.sh ../data/functional-protein-exon2-dataset.csv "$HISTONE_NAME" protein-exon2-histone-feature
sh histone_processing.sh ../data/protein-exon2-negative-control-dataset-curated.csv "$HISTONE_NAME" protein-exon2-NC-histone-feature

sh histone_processing.sh ../data/functional-protein-exon3-dataset.csv "$HISTONE_NAME" protein-exon3-histone-feature
sh histone_processing.sh ../data/protein-exon3-negative-control-dataset-curated.csv "$HISTONE_NAME" protein-exon3-NC-histone-feature

sh join_datasets.sh ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_protein-exon2-histone-feature.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_protein-exon3-histone-feature.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_protein_positives_matrix
sh join_datasets.sh ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_protein-exon2-NC-histone-feature.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_protein-exon3-NC-histone-feature.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_protein_negatives_matrix
sh join_datasets.sh ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_protein_positives_matrix.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_protein_negatives_matrix.csv ../data/histone_feature/"$HISTONE_NAME"/"$HISTONE_NAME"_protein_matrix

