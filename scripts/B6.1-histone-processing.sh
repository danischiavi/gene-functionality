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
output_path="data/epigenetic/histone_feature/${HISTONE_NAME}"
marks_path="data/raw/epigenetic_data/histone_marks/${HISTONE_NAME}"

# Create a directory if not present already
if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi

# Remove header from regions file and sort regions
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"


# Heavy processing ahead. Check for previously created files and skip if present.
if [ ! -f "${marks_path}/${HISTONE_NAME}_sorted_histone_marks.bed" ]; then
    for dir in "$marks_path"/*; do
    	if [ -d "$dir" ]; then 
			for file in "$dir"/*.bed.gz; do
				if [ -f "$file" ]; then
        			filename=$(basename "$file" .bed)  # Extract filename without .bed extension
        			zcat "$file" >> "${marks_path}/temp_histone_marks.bed"
					#awk -v OFS="\t" '{print $0}' "$file" >> "$marks_path"/temp_histone_marks.bed
				fi
			done
		fi
    done

    # Select columns of interest from marks file, i.e. Chromosome, Start, End, Mark Signal Value, p-value, and q-value.
    # Then sort by chromosome, start, and end.
    awk -F '\t' '{
             print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9
         }' OFS=, "${marks_path}/temp_histone_marks.bed" | \
    sort -k1,1 -k2,2n -k3,3n > "${marks_path}/${HISTONE_NAME}_sorted_histone_marks.bed"
fi

# Perform interval intersection using BEDTools
if [[ ! -f "${OUTPUT_NAME}-overlaps.bed" ]]; then
	bedtools intersect -a "${marks_path}/${HISTONE_NAME}_sorted_histone_marks.bed" -b "$sorted_regions_file_bed" -wo > "${OUTPUT_NAME}-overlaps.bed"
fi

# AvgSignal = sum(wi * ei)
# Calculate the weighted sum of enrichment, grouped by Chromosome, StartExon, EndExon
awk -F '\t' '
    BEGIN { OFS="\t" }
    {
        # Columns from -wo output:
        # 1-6: Fields from file A ( histone marks: chr, start, end, signal, pval, qval)
        # 7-11: Fields from file B ( regions: chr, start, end, ID, Functional)
        # 12: Length of overlap

        # Define the key for the region (Chromosome, StartRegion, EndRegion from file B)
        region_chr = $7;
        region_start = $8;
        region_end = $9;
        key = region_chr FS region_start FS region_end;

        # Get values for calculation
        signal = $4;         # Signal value from histone mark (col 4)
        overlap_len = $12;   # Length of the overlap (col 12)
        region_len = region_end - region_start; # Length of the region from file B

        # Avoid division by zero if region length is 0 (should not happen with valid BED)
        if (region_len > 0) {
            # Calculate weighted signal component for this overlap (signal * overlap_fraction)
            weighted_signal_component = signal * (overlap_len / region_len);
            sum_weighted_signal[key] += weighted_signal_component;

            # Calculate scaled signal for this overlap (signal * overlap / region_length)
            # Note: This is the same as weighted_signal_component in this context
            scaled_signal = weighted_signal_component;

            # Update maximum scaled signal if current scaled signal is higher
            if (!(key in max_scaled_signal) || scaled_signal > max_scaled_signal[key]) {
                max_scaled_signal[key] = scaled_signal;
            }
        } else {
            # Handle potential zero-length regions if necessary
            print "Warning: Zero length region encountered:", $7, $8, $9 > "/dev/stderr";
        }
    }
    END {
        # Print the results for each region
        # Output format: Chr, StartRegion, EndRegion, AvgSignal, MaxSignal, MaxScaledSignal
        for (k in sum_weighted_signal) {
            
            # Get max values, defaulting to 0 if not found (e.g., if only zero-length regions overlapped)
            max_scaled_sig = (k in max_scaled_signal) ? max_scaled_signal[k] : 0;

            print k, max_scaled_sig;
        }
    }
' "${OUTPUT_NAME}-overlaps.bed" | sort -k1,1 -k2,2n -k3,3n > "${OUTPUT_NAME}.bed"


awk -F '\t' 'NR==FNR { seen[$1,$2,$3]=$4; next } 
             { if (($1,$2,$3) in seen) print $0"\t"seen[$1,$2,$3]; else print $0"\t0" }' \
             "${OUTPUT_NAME}.bed" "$sorted_regions_file_bed" > "$augmented_regions_bed"

# Add header row
echo "${HISTONE_NAME}_MaxScaledSignal" > "${OUTPUT_NAME}.csv"

# Convert BEDTools output to desired CSV format
awk -F '\t' '{print $6}' OFS=, "$augmented_regions_bed" >> "${OUTPUT_NAME}.csv"


# Remove the temporary files
# rm "$temp_histone_marks_bed"
# rm "$sorted_histone_marks_bed"
rm "$sorted_regions_file_bed"
#rm "$overlaps_bed"
rm "$augmented_regions_bed"
rm "$OUTPUT_NAME".bed
