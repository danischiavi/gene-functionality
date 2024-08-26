#!/bin/bash
#
# Dependencies: 
# - access_py.py: https://github.com/bkb3/openen/tree/master
#
# Script Name: accessibility_processing.sh
#
# Author: Estefania Rojas
#
# Description: This script calculates accessibility scores for sequences of functional long and short-ncRNA 
# databases and protein-coding-RNA databases.
#
# Input: 1. Csv file containing the sequeces dataset to be analysed. Ex. data/functional-lncrna-exon1-dataset.csv
#        2. Name of the output csv file to save accessibility scores for each sequence in input file 1.
#           Ex. lncrna-exon1-accessibility-feature
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
output_path=../data/accessibility_feature

mkdir -p "$output_path"

# Print output header
echo "accessibility" > "$output_path"/"$OUTPUT_NAME".csv

# Remove header from regions file and get just the sequences column
awk -F, 'NR > 1 {print $6}' OFS="\t", "$REGIONS_FILE" | while read -r line
do
    python3 B0.4.1-access_py.py -s "$line" >> "$output_path"/"$OUTPUT_NAME".csv
done
