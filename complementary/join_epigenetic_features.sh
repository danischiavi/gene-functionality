#!/bin/bash

set -uex

# Check for correct usage
if [ $# -ne 3 ]; then
    echo "Usage: $0 initial_matrix.csv list_of_files.txt output_matrix_name"
    exit 1
fi

# Store input files with descriptive names
first_file=$1
list_of_files=$2
output=$3

# Extract common columns
cut -d $'\t' -f 5-6 "$first_file" > temp2_file.csv

for file in $(cat "$list_of_files"); do
    cut -d $'\t' -f6 $file > temp_file.csv
    paste -d $'\t' temp2_file.csv temp_file.csv > temp3_file.csv
    mv temp3_file.csv temp2_file.csv     # Update temp2 file
done


echo Functionality$'\t'chrm_acc$'\t'H3k27ac$'\t'H3k27me3$'\t'H3k36me3$'\t'H3k4me1$'\t'H3k4me3$'\t'H3k9me3$'\t'methylome_feature > "$output"

grep -v "Functional" temp2_file.csv | awk -F'\t' '{print $0 }' >> "$output"



# Remove temps
rm -rf temp_file.csv
rm -rf temp2_file.csv
rm -rf temp2_file.csv
