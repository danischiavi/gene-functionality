#!/bin/bash
#
# Dependencies: 
# - bedtools:
#       https://bedtools.readthedocs.io/en/latest/content/installation.html
# - tabix:
#		https://www.htslib.org/doc/tabix.html
# - bigWigToBedGraph:
#		https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
#
# Script Name: GERP_processing.sh
#
# Author: Estefania Rojas
#
# Description: This script calculates the mean and max gerp score for sequences of functional long and short-ncRNA 
# databases and protein-coding-RNA databases. This is done using 'bedtools map' command.
#
# Input: 1. Csv file containing the sequeces dataset to be analysed. Ex. data/functional-lncrna-exon1-dataset.csv
#        2. bw file containing gerp scores
#        3. Name of the output csv file to save mean and max features for each sequence in input file 1.
#           Ex. lncrna-exon1-histone-feature
#
############################################################################################################################

# Bash strict mode. -e exit on error, -u error on unset variables, -x print every command executed.
set -ue

# Get the date of execution. Format example: Tuesday January 17, 2024
DATE=$(date "+%A %B %d, %Y")

# Print the date to console
echo ${DATE}
##############################################################

# Check for correct usage
if [ $# -ne 4 ]; then
    echo "Usage: $0 regions_to_extract.csv gerp_scores_file.bedgraph prefix"
    exit 1
fi

# Store input files with descriptive names
REGIONS_FILE=$1
GERP_FILE=$2
PREFIX=$3

output_directory=data/conservation
mkdir -p "$output_directory"

output_file="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')gerp${PREFIX}-feature.csv"


# Temporary files
sorted_regions_file_bed=$(mktemp)

# Remove header from regions file, sort regions and convert to bed format
awk -F, 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1"\t"$2}' OFS="\t", "$REGIONS_FILE" | \
sort -k1,1 -k2,2n -k3,3n > "$sorted_regions_file_bed"

# Check if sorted bedgraph doesnt exist
if [ ! -f "${GERP_FILE%.bw}"_sorted.bedgraph ]; then
	# Convert file from bigWig to bedGraph using bigWigToBedGraph (UCSC Kent Tools)
	bigWigToBedGraph "$GERP_FILE" "${GERP_FILE%.bw}".bedgraph

	# Order bedgraph file, use awk to add text to the first column of each row
	sort -k1,1 -k2,2n "${GERP_FILE%.bw}".bedgraph | tr ' ' '\t' | \
	awk -v prefix="chr" 'BEGIN {FS=OFS="\t"} {$1=prefix $1; print $0}' > "${GERP_FILE%.bw}"_sorted.bedgraph
fi

# Compute max gerp score
bedtools map -a "$sorted_regions_file_bed" -b "${GERP_FILE%.bw}"_sorted.bedgraph -c 4 -o max -sorted | \
awk -F '\t' '{print $6"\t"$7}' > "$output_file"_augmented

# Add header row to output file
echo GERP_"$PREFIX"_max > "$output_file"

# Convert BED output to desired CSV format
sed 's/\t/,/g' "$output_file"_augmented >> "$output_file"

rm "$sorted_regions_file_bed"
rm "$output_file"_augmented
