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
#temp_chrm_acc_marks_bed=$(mktemp)  # Temporary file for BED conversion
#sorted_chrm_acc_marks_bed=$(mktemp)
sorted_regions_file_bed=$(mktemp)
#overlaps_bed=$(mktemp)
augmented_regions_bed=$(mktemp)

# Variables
data_path=data/epigenetic/chrm_acc_feature
marks_path=data/raw/epigenetic_data/chromatin_accessibility

# Create a directory if not present already
if [ ! -d "$data_path" ]; then
    mkdir "$data_path"
fi


# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"


# Heavy processing ahead. Check for previously created files and skip if present.
if [ ! -f "${marks_path}/chrm_acc_sorted_marks.bed" ]; then
    for dir in "$marks_path"/*; do
		if [ -d "$dir" ]; then 
			for file in "$dir"/*.bed.gz; do
				if [ -f "$file" ]; then
					filename=$(basename "$file" .bed)  # Extract filename without .bed extension
        			zcat "$file" >> "${marks_path}/temp_histone_marks.bed"
				fi
			done
		fi
	done     

    awk -F '\t' '{
             print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9
         }' OFS=, "${marks_path}/temp_chrm_acc_marks.bed" | \
    sort -k1,1 -k2,2n -k3,3n > "${marks_path}/chrm_acc_sorted_marks.bed"
fi


# Perform interval intersection using BEDTools
if [[ ! -f "${OUTPUT_NAME}-overlaps.bed" ]]; then
bedtools intersect -a "${marks_path}/chrm_acc_sorted_marks.bed" -b "$sorted_regions_file_bed" -wo > "${OUTPUT_NAME}-overlaps.bed"
fi

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
' "${OUTPUT_NAME}-overlaps.bed" | sort -k1,1 -k2,2n -k3,3n > "${OUTPUT_NAME}.bed"


awk -F '\t' 'NR==FNR { seen[$1,$2,$3]=$4; next } 
             { if (($1,$2,$3) in seen) print $0"\t"seen[$1,$2,$3]; else print $0"\t0" }' \
             "${OUTPUT_NAME}.bed" "$sorted_regions_file_bed" > "$augmented_regions_bed"

# Add header row
echo "chrm_acc_AvgSignal" > "${OUTPUT_NAME}.csv"

# Convert BEDTools output to desired CSV format
awk -F '\t' '{print $6}' OFS=, "$augmented_regions_bed" >> "${OUTPUT_NAME}.csv"


# Remove the temporary files
# rm "$temp_chrm_acc_marks_bed"
# rm "$sorted_chrm_acc_marks_bed"
rm "$sorted_regions_file_bed"
#rm "$overlaps_bed"
rm "$augmented_regions_bed"
rm "$OUTPUT_NAME".bed
