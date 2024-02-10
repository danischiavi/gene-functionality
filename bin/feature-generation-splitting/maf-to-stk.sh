#!/bin/bash
#
# Convert MAF (multiple alignment format) to Stockholm format: 
#       MAF blocks for each species are concatenated in the final stk file
#       Sequences which don't match human sequence length are filter out ( $len_38 can be modify for other reference seq)
#
# INPUT: 
#       $1 MAF file for conversion
#       $2 directory to store output stk format 
#
# Author: Daniela Schiavinato - February 2024 
#
############################################################################################################################

## Files and directory

input_file=$1
stk_directory=$2

id="$(basename "${input_file%.*}" | sed 's/maf//')"
output_file="${stk_directory}/$(basename "${input_file%.*}" | sed 's/maf//').stk"


## Array with species Ids present in the MAF alignment 

declare -a "species_id=()"
mapfile -t species_id < <(grep -v -E "#\|score|^$" "$input_file" | awk '{print $2}' | awk '!seen[$0]++')    


## Find the longest Id to align sequences in stk file relatively to it

longest_id=$( 
    echo "${species_id[@]}" | \
    tr ' ' '\n' | \
    awk '{ if (length > max) { max = length; longest = $0 } } END { print length(longest) }' 
)


## Extract sequences from MAF file following the IDs and append them on a new file

for species_id in "${species_id[@]}"; do

    sequence=$( awk -v id="$species_id" '$2 ~ id {print $NF}' "$input_file" | awk '{printf "%s", $0}' ) 
    printf "%-*s\t%s\n" "$longest_id" "$species_id"  "$sequence" >> "$output_file".tmp    

done


## stk file - print header; remove partial aligments (seqs which don't match with human seqs length)

len_hg38=$(awk -v partial_id="hg38" '{ if ($1 ~ "^" partial_id) { print length($2) } }' "$output_file".tmp)

echo -e "# STOCKHOLM 1.0 \n#=GF ID\t$id" > "$output_file"
awk -v len="$len_hg38" 'length($2) == len' "$output_file".tmp > "$output_file"
echo "//" >> "$output_file" 

rm -rf "$output_file".tmp

############################################################################################################################
