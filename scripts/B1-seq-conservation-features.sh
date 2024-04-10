#!/bin/bash
#
# Script Name: seq-conservation-features.sh
#
# Description: This script calculates the PhyloP (based on 241way-alignment of Zoonomia project), PhastCons and GERP (based on 100way vertebrate alignment) for for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences.
#
# Input: $1 is the dataset 
#       
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.
#
############################################################################################################################

#### Files and directories #### 
initial_data=$1
bigWigSummary_exe=$2
zoonomia_phylo_bw=$4

output_directory=data/conservation
mkdir -p "$output_directory"

output_file="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')conservation.csv"                  


#### Extract conservation features for each set of chromosome coordinates #### 
if [ ! -s "$output_file" ]; then

    echo "phyloP_max" > "$output_file"

    tail -n +2 "$initial_data"  | while IFS=',' read -r _ _ chr start end _; do

        #### PhyloP (pp) values: 241-way mammalian alignment ####
        echo "zoonomia_max=$(   $bigWigSummary_exe -type=max $zoonomia_phylo_bw $chr $start $end 1 2>&1 )" >> errors.log
        
        zoonomia_max=$(   $bigWigSummary_exe -type=max "$zoonomia_phylo_bw" "$chr" "$start" "$end" 1 2>&1 )

        ## Convert any missing data to NA
        missing_value() {

            local var="$1"

            test_count=$(echo "$var" | wc -w)
            
            if [ "$test_count" -ne 1 ]; then var="NA"; fi
        
            echo "$var"
        }

        zoonomia_max=$(missing_value "$zoonomia_max")

        ## Values to output file 
        echo "$zoonomia_max" >> "$output_file"

    done

fi



