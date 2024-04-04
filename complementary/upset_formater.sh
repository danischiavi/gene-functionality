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
cut -d $'\t' -f 1-3 "$first_file" > temp2_file.csv
# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $0}' OFS="\t", temp2_file.csv | \
sort -k1,1 -k2,2n -k3,3n > sorted_regions_file.bed

mv sorted_regions_file.bed temp2_file.csv

for file in $(cat "$list_of_files"); do
    if [ "$first_file" == "$file" ]; then
        # Fill column of ones when intersection is with itself
        awk -F '\t' '{print $0"\t1"}' OFS="\t", temp2_file.csv > sorted_regions_file.bed
        mv sorted_regions_file.bed temp2_file.csv
    fi
    if [ ! "$first_file" == "$file" ]; then
        # Remove header from database file and sort regions
        awk -F, 'NR > 1 {print $0}' OFS="\t", $file | \
        sort -k1,1 -k2,2n -k3,3n > sorted_db.bed
        
        bedtools intersect -f 0.9 -F 0.9 -e -a temp2_file.csv -b sorted_db.bed -wa -u | sort -k1,1 -k2,2n -k3,3n > found.bed

        awk -F '\t' 'NR==FNR { seen[$1,$2,$3]=1; next } 
             { if (($1,$2,$3) in seen) print $0"\t"seen[$1,$2,$3]; else print $0"\t0" }' \
             found.bed temp2_file.csv > augmented.bed

        mv augmented.bed temp2_file.csv
    fi

done


echo Chromosome$'\t'Start$'\t'End$'\t'GencodeV44$'\t'GencodeV45$'\t'HGNC$'\t'NCBI$'\t'RNACentral$'\t'UniProt > "$output"

awk -F'\t' '{print $0 }' temp2_file.csv >> "$output"



# Remove temps
rm -rf sorted_regions_file.bed
rm -rf temp2_file.csv
rm -rf found.bed
rm -rf augmented.bed
