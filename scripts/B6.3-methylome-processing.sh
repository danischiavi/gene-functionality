#!/bin/bash

set -uex

# Check for correct usage
if [ $# -ne 2 ]; then
    echo "Usage: $0 regions_to_extract.csv outname"
    exit 1
fi

# Store input files with descriptive names
REGIONS_FILE=$1
OUTPUT_NAME=$2

# Temporary files
sorted_regions_file_bed=$(mktemp)
augmented_regions_bed=$(mktemp)

# Variables
data_path=data/raw/epigenetic_data/methylome_feature
beds_path=data/raw/epigenetic_data/methylome/methylome_beds

# Create a directory if not present already
if [ ! -d "$data_path"/sorted_beds/ ]; then
    mkdir -p "$data_path"/sorted_beds/
fi

# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"

for file in "$beds_path"/*.bed; do
    filename=$(basename "$file" .bed)
    if [ ! -f "$data_path"/sorted_beds/"$filename"_sorted.bed ]; then
        #Get columns chrom, start, end, methylation level, and percentage of methylation
        awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$10"\t"$11}' "$file" | \
        sort -k1,1 -k2,2n -k3,3n > "$data_path"/sorted_beds/"$filename"_sorted.bed
    fi

    # Perform interval intersection using BEDTools
    bedtools intersect -a "$data_path"/sorted_beds/"$filename"_sorted.bed -b "$sorted_regions_file_bed" -wa -wb >> overlaps_"$OUTPUT_NAME".bed
done

# Group hits by chromosome, start, and end of reads to obtain average methilation percentage
awk -F '\t' '
    BEGIN { OFS="\t" }
    {
        key = $1 FS $7 FS $8;
        sum[key]+=$5;
        count[key]++
    }
    END {
        for (k in sum) {
            print k, sum[k]/count[k];
        }
    }
' overlaps_"$OUTPUT_NAME".bed | sort -k1,1 -k2,2n -k3,3n > "$OUTPUT_NAME".bed


awk -F '\t' 'NR==FNR { seen[$1,$2,$3]=$4; next } 
             { if (($1,$2,$3) in seen) print $0"\t"seen[$1,$2,$3]; else print $0"\t0" }' \
             "$OUTPUT_NAME".bed "$sorted_regions_file_bed" > "$augmented_regions_bed"


# Add header row
echo Chromosome'\t'Start'\t'End'\t'ID'\t'Functional'\t'methylAvgPercentage > "$data_path"/"$OUTPUT_NAME".csv

# Convert BEDTools output to desired CSV format
awk -F '\t' '{print $0}' OFS=, "$augmented_regions_bed" >> "$data_path"/"$OUTPUT_NAME".csv


# Remove the temporary files
rm "$sorted_regions_file_bed"
rm "$augmented_regions_bed"
rm "$OUTPUT_NAME".bed

