#!/bin/bash
#
# Dependencies: 
# - markovProperties.pl: Used to compute the base 2 log of the dinucleotide frequences in a given sequence.
#
# Script Name: dinucleotide_frequencies.sh
#
# Author: Estefania Rojas
#
# Description: This script calculates dinucleotide frequencies for sequences of functional long and short-ncRNA 
#              databases and protein-coding-RNA databases for all possible dincleotides. The result is savet at
#              the directory specified by $output_path variable.
#
# Input: 1. Csv file containing the sequeces dataset to be analysed. Ex. data/functional-lncrna-exon1-dataset.csv
#        2. Name of the output csv file to save dinucleotide frequencies for each sequence in input file 1.
#           Ex. lncrna-exon1-dinucleotide-features
#
############################################################################################################################

# Bash strict mode. -e report all errors, -u exit on error, -x print every command executed.
set -uex

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
output_path=../data/datasets/dinucleotide_feature

mkdir -p "$output_path"

# Remove header from regions file and get just the sequences column
awk -F, 'NR > 1 {print $6}' OFS="\t", "$REGIONS_FILE" > sequences.fa

# Print output header
echo AA,AC,AG,AT,CA,CC,CpG,CT,GA,GC,GG,GT,TA,TC,TG,TT > dinucleotide_freqs.csv

#cat sequences.fa | while read -r line
#do
#    python compute_dinucleotide_log2freq.py "$line" >> dinucleotide_freqs.csv
#done
cat sequences.fa | while read -r line
do
    echo "$line" > sequence.fa
    ./markovProperties.pl -k 2 -i sequence.fa >> dinucleotide_freqs.csv
    #echo "" >> dinucleotide_freqs.csv
done

awk -F, '{print $0}' OFS="," dinucleotide_freqs.csv > "$output_path"/"$OUTPUT_NAME".csv

rm dinucleotide_freqs.csv
rm sequences.fa
rm sequence.fa
