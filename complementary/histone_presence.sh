#!/bin/bash

set -uex

# Check for correct usage
if [ $# -ne 3 ]; then
    echo "Usage: $0 regions_to_compare.csv threshold outname"
    exit 1
fi

# Store input files with descriptive names
regions_file=$1
threshold=$2
outname=$3

# Temporary files
temp_histone_marks_bed=$(mktemp)  # Temporary file for BED conversion
sorted_histone_marks_bed=$(mktemp)
sorted_regions_file_bed=$(mktemp)
overlaps_bed=$(mktemp)


# Combine all bed files with H3 histone marks into one.
cat ../data/histone_marks/H3k36me3/*.bed > "$temp_histone_marks_bed"


# Calculate -log(threshold) for filtering by p and q values later
log_threshold=$(awk -v thres="$threshold" 'BEGIN {print -log(thres)}')

#cut -f1,2,3,7 "$temp_histone_marks_bed" | \

# Filtering based on p and q value (-log10 format), select only marks woth values less than 0.05
awk -F '\t' -v thres="$log_threshold" '
         ($8 != -1 && $8 > thres || $8 == -1 ) &&
         ($9 != -1 && $9 > thres || $9 == -1 )  {
         print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9
        }' OFS=, "$temp_histone_marks_bed" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_histone_marks_bed"

# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$regions_file" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"

# Perform interval intersection using BEDTools
bedtools intersect -a "$sorted_regions_file_bed" -b "$sorted_histone_marks_bed" -wa -u > overlaps.bed

awk -F '\t' 'NR==FNR { seen[$1,$2,$3]="yes"; next } 
             { if (($1,$2,$3) in seen) print $0,"yes"; else print $0,"no" }' \
             overlaps.bed "$sorted_regions_file_bed" > augmented_regions.bed

# Add header row
echo Chromosome,Start,End,ID,Functional,HistonePresence > "$outname".csv
cat augmented_regions.bed >> "$outname".csv
