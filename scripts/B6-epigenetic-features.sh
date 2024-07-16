#!/bin/bash

initial_data=$1

output_directory=data/epigenetic
mkdir -p "$output_directory"

file_name=${output_directory}/$(basename "${initial_data%.*}" | sed 's/-dataset//')

# Output file
output_file_histone="${file_name}-histones"
output_file="$file_name"-epigenetic.csv 


#### Run scripts #### 
# Run histone signatures script
histone_name=("H3K36me3" "H3K27ac" "H3K79me2")

if [ ! -s "$output_file_histone" ]; then	

	for histone in "${histone_name[@]}"; do
	
		./scripts/B6.1-histone-processing.sh "$initial_data" "$histone" "${file_name}-${histone}" 

	done

	paste -d',' "${file_name}-H3K36me3" "${file_name}-H3K27ac" "${file_name}-H3K79me2" > "$output_file_histone"

	# rm -rf "${file_name}-H3K36me3" "${file_name}-H3K27ac" "${file_name}-H3k79me2"

fi

# Run chromatin access 
if [ ! -s "${file_name}-chrm-acc" ]; then

	./scripts/B6.2-chrm-acc-processing.sh "$initial_data" "${file_name}-chrm-acc" 

fi


# Run methylome signature 
if [ ! -s "${file_name}-methylome" ]; then

	./scripts/B6.3-methylome-processing.sh "$initial_data" "${file_name}-methylome" 

fi



## Join output files for better organization into one
if [ ! -s "$output_file_specific" ]; then

    paste -d',' "$output_file_histone" "${file_name}-chrm-acc" "${file_name}-methylome" > "$output_file"

fi