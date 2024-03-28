#!/bin/bash

set -uex


# LONG NCRNA
## Process histone marks for exon1
sh chrm_acc_processing.sh ../data/functional-lncrna-exon1-dataset.csv lncrna-exon1-histone-feature
sh chrm_acc_processing.sh ../data/lncrna-exon1-negative-control-dataset-curated.csv lncrna-exon1-NC-histone-feature
## Process histone marks for exon2
sh chrm_acc_processing.sh ../data/functional-lncrna-exon2-dataset.csv lncrna-exon2-histone-feature
sh chrm_acc_processing.sh ../data/lncrna-exon2-negative-control-dataset-curated.csv lncrna-exon2-NC-histone-feature
## Join datasets into one for R script
sh join_datasets.sh ../data/chrm_acc_feature/chrm_acc_lncrna-exon1-histone-feature.csv ../data/chrm_acc_feature/chrm_acc_lncrna-exon2-histone-feature.csv ../data/chrm_acc_feature/chrm_acc_lncrna_positives_matrix
sh join_datasets.sh ../data/chrm_acc_feature/chrm_acc_lncrna-exon1-NC-histone-feature.csv ../data/chrm_acc_feature/chrm_acc_lncrna-exon2-NC-histone-feature.csv ../data/chrm_acc_feature/chrm_acc_lncrna_negatives_matrix
sh join_datasets.sh ../data/chrm_acc_feature/chrm_acc_lncrna_positives_matrix.csv ../data/chrm_acc_feature/chrm_acc_lncrna_negatives_matrix.csv ../data/chrm_acc_feature/chrm_acc_lncrna_matrix

#SHORT NCRNA
sh chrm_acc_processing.sh \
    ../data/functional-short-ncrna-dataset.csv \
    \
    short-ncrna-histone-feature
sh chrm_acc_processing.sh \
    ../data/short-ncrna-negative-control-dataset-curated.csv \
    \
    short-ncrna-NC-histone-feature
sh join_datasets.sh \
    ../data/chrm_acc_feature/chrm_acc_short-ncrna-histone-feature.csv \
    ../data/chrm_acc_feature/chrm_acc_short-ncrna-NC-histone-feature.csv \
    ../data/chrm_acc_feature/chrm_acc_short_ncrna_matrix


# PROTEIN CODING
sh chrm_acc_processing.sh ../data/functional-protein-exon2-dataset.csv protein-exon2-histone-feature
sh chrm_acc_processing.sh ../data/protein-exon2-negative-control-dataset-curated.csv protein-exon2-NC-histone-feature

sh chrm_acc_processing.sh ../data/functional-protein-exon3-dataset.csv protein-exon3-histone-feature
sh chrm_acc_processing.sh ../data/protein-exon3-negative-control-dataset-curated.csv protein-exon3-NC-histone-feature

sh join_datasets.sh ../data/chrm_acc_feature/chrm_acc_protein-exon2-histone-feature.csv ../data/chrm_acc_feature/chrm_acc_protein-exon3-histone-feature.csv ../data/chrm_acc_feature/chrm_acc_protein_positives_matrix
sh join_datasets.sh ../data/chrm_acc_feature/chrm_acc_protein-exon2-NC-histone-feature.csv ../data/chrm_acc_feature/chrm_acc_protein-exon3-NC-histone-feature.csv ../data/chrm_acc_feature/chrm_acc_protein_negatives_matrix
sh join_datasets.sh ../data/chrm_acc_feature/chrm_acc_protein_positives_matrix.csv ../data/chrm_acc_feature/chrm_acc_protein_negatives_matrix.csv ../data/chrm_acc_feature/chrm_acc_protein_matrix

