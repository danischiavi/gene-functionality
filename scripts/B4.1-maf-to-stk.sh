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
mapfile -t species_id < <(grep -v -E "#|score=[[:alnum:]]|^$" "$input_file" | awk '{print $2}' | awk '!seen[$0]++')    

species_id_count=${#species_id[@]}

if [ "$species_id_count" -gt 1 ]; then 

	##### Find the longest Id to align sequences in stk file relatively to it
	longest_id=$( 
    	echo "${species_id[@]}" | \
    	tr ' ' '\n' | \
    	awk '{ if (length > max) { max = length; longest = $0 } } END { print length(longest) }' 
	)


	#### Extract sequences from MAF file following the IDs and append them on a new file
	for species_id in "${species_id[@]}"; do

    	sequence=$( awk -v id="$species_id" '$2 ~ id {print $NF}' "$input_file" | awk '{printf "%s", $0}' ) 
    	printf "%-*s\t%s\n" "$longest_id" "$species_id"  "$sequence" >> "${output_file}.tmp"    

	done


	## Remove partial alignments < 50%: determinated by length less than 50% of human's not including '-' and 'N'
	# Pad included partial alignments with 'N' to match human's len for further analysis with tools (RNAcode, R-scape and RNAalifold) 

	len_hg38=$(awk -v partial_id="hg38" '{ if ($1 ~ "^" partial_id) { print length($2) } }' "${output_file}.tmp")

	awk -v len="$len_hg38" -v longest_id="$longest_id" '
    	function custom_length(str) {
        	gsub("-", "", str);                                         
        	gsub("N", "", str) 
        	return length(str); 
    	} 

    	function pad_sequence(seq, target_len, longest_id) {
        	sequence_length = length(seq)
        	padding_length = target_len - sequence_length
        	if (padding_length > 0) {
            	padding = sprintf("%*s", padding_length, "")
            	gsub(/ /, "N", padding)
            	return seq padding
        	} else {
            	return seq
        	}
    	}

    	{ 
        	if (custom_length($NF) >= len * 0.5 && length($NF) <= len) {
            	sequence = pad_sequence($NF, len, longest_id)
            	printf "%-*s\t%s\n", longest_id, $1, sequence
        	}   
    	}' "${output_file}.tmp" > "${output_file}.tmp2"


	## Format stk file
	echo -e "# STOCKHOLM 1.0 \n#=GF ID\t$id" > "$output_file"
	awk '!seen[$2]++' "${output_file}.tmp2" >> "$output_file"         # Remove duplicated sequences
	echo "//" >> "$output_file" 

	rm -rf "${output_file}.tmp" "${output_file}.tmp2"
	
	count_lines=$( wc -l < "$output_file" )

	if [ "$count_lines" -le 3 ]; then 
		echo "maf file NA" > "$output_file"
	fi

else 

	echo "maf file NA" > "$output_file" 

fi

############################################################################################################################
