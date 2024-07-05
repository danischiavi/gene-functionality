#!/bin/bash

set -uex

# Check for correct usage
if [ $# -ne 3 ]; then
    echo "Usage: $0 regions_to_extract.csv histone_name outname"
    exit 1
fi

# Store input files with descriptive names
REGIONS_FILE=$1
HISTONE_NAME=$2
OUTPUT_NAME=$3

# Temporary files
#temp_histone_marks_bed=$(mktemp)  # Temporary file for BED conversion
#sorted_histone_marks_bed=$(mktemp)
sorted_regions_file_bed=$(mktemp)
#overlaps_bed=$(mktemp)
augmented_regions_bed=$(mktemp)

# Variables
marks_path=data/raw/histone_marks/"$HISTONE_NAME"


# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"


# Heavy processing ahead. Check for previously created files and skip if present.
if [ ! -f "$data_path"/"$HISTONE_NAME"_sorted_histone_marks.bed ]; then
    for file in "$marks_path"/*.bed; do
        filename=$(basename "$file" .bed)  # Extract filename without .bed extension
        awk -v OFS="\t" '{print $0}' "$file" >> "$marks_path"/temp_histone_marks.bed
    done

    # Filtering based on p and q value (-log10 format), select only marks with values less than 0.05
    awk -F '\t' '{
             print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9
         }' OFS=, "$marks_path"/temp_histone_marks.bed | \
    sort -k1,1 -k2,2n -k3,3n > "$data_path"/"$HISTONE_NAME"_sorted_histone_marks.bed
fi


# Perform interval intersection using BEDTools
bedtools intersect -a "$data_path"/"$HISTONE_NAME"_sorted_histone_marks.bed -b "$sorted_regions_file_bed" -wo > overlaps_"$HISTONE_NAME"_"$OUTPUT_NAME".bed

# AvgSignal = sum(wi * ei)
# Calculate the weighted sum of enrichment, grouped by Chromosome, StartExon, EndExon
awk -F '\t' '
    BEGIN { OFS="\t" }
    {
        key = $1 FS $8 FS $9;
        sum[key] += $4 * ($12 / ($9 - $8));
    }
    END {
        for (k in sum) {
            print k, sum[k];
        }
    }
' overlaps_"$HISTONE_NAME"_"$OUTPUT_NAME".bed | sort -k1,1 -k2,2n -k3,3n > "$OUTPUT_NAME".bed


awk -F '\t' 'NR==FNR { seen[$1,$2,$3]=$4; next } 
             { if (($1,$2,$3) in seen) print $0"\t"seen[$1,$2,$3]; else print $0"\t0" }' \
             "$OUTPUT_NAME".bed "$sorted_regions_file_bed" > "$augmented_regions_bed"

# Add header row
echo Chromosome'\t'Start'\t'End'\t'ID'\t'Functional'\t'AvgSignal > "$data_path"/"$HISTONE_NAME"_"$OUTPUT_NAME".csv

# Convert BEDTools output to desired CSV format
awk -F '\t' '{print $0}' OFS=, "$augmented_regions_bed" >> "$data_path"/"$HISTONE_NAME"_"$OUTPUT_NAME".csv


# Remove the temporary files
# rm "$temp_histone_marks_bed"
# rm "$sorted_histone_marks_bed"
rm "$sorted_regions_file_bed"
#rm "$overlaps_bed"
rm "$augmented_regions_bed"
rm "$OUTPUT_NAME".bed
