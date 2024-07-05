#!/bin/bash

set -uex

# Check for correct usage
if [ $# -ne 3 ]; then
    echo "Usage: $0 positive_dataset negative_dataset output_matrix_name"
    exit 1
fi

# Store input files with descriptive names
dataset1=$1
dataset2=$2
output=$3

# Temporary files
temp_positive=$(mktemp)  
temp_negative=$(mktemp)

# Remove header
grep -v "Start" "$dataset1" | awk -F'\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5"\t" $6 }' > "$temp_positive"
grep -v "Start" "$dataset2" | awk -F'\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5"\t" $6 }' > "$temp_negative"

# Add header row
echo Chromosome'\t'Start'\t'End'\t'ID'\t'Functional'\t'AvgSignal > "$output".csv

# join datasets
cat "$temp_positive" "$temp_negative" | \
sort -k2,2 -k3,3n -k4,4n >> "$output".csv

rm "$temp_positive"
rm "$temp_negative"
