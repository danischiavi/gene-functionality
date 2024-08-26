#!/bin/bash
#
# Dependencies: 
# - CPC2.py: https://github.com/gao-lab/CPC2_standalone
#
# Script Name: fickett_score_processing.sh
#
# Author: Estefania Rojas
#
# Description: This script calculates fickett score for sequences of functional long and short-ncRNA 
# databases and protein-coding-RNA databases.
#
# Input: 1. Csv file containing the sequeces dataset to be analysed. Ex. data/functional-lncrna-exon1-dataset.csv
#        2. Name of the output csv file to save fickett score for each sequence in input file 1.
#           Ex. lncrna-exon1-fickett-features
#
############################################################################################################################

# Bash strict mode. -e report all errors, -u exit on error, -x print every command executed.
set -ue

# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}
##############################################################

# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 sequences_dataset.csv outname"
    exit 1
fi

# Store input files with descriptive names
REGIONS_FILE=$1
OUTPUT_NAME=$2

# Define path to store output
output_path=../data/fickett_feature

mkdir -p "$output_path"

# Remove header from regions file and get just the sequences column simulate fasta format
awk -F, 'NR > 1 {print ">"$1"\n"$6}' OFS="\t", "$REGIONS_FILE" > sequences.fa

# Compute fickett score
python3 bin/CPC2_standalone-1.0.1/bin/CPC2.py -i sequences.fa -o fickett_scores

# Print output header
#echo Fickett_score > "$output_path"/"$OUTPUT_NAME".csv

# keep just Ficket_score column from CPC2.py output
awk -F '\t' '{print $4}' OFS="," fickett_scores.txt > "$output_path"/"$OUTPUT_NAME".csv


#rm -rf functional.csv
#rm -rf result.csv
rm -rf fickett_scores.txt
#rm -rf sequence.fa
rm -rf sequences.fa
