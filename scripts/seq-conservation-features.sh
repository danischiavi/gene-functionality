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
phast_bw=$3
zoonomia_phylo_bw=$4
gerp_bw=$5

output_directory=data/conservation
mkdir -p "$output_directory"

output_file="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')conservation.csv"                  


#### Extract conservation features for each set of chromosome coordinates #### 
if [ ! -s "$output_file" ]; then

    echo "phyloP_mean,phyloP_max, phastCons_mean,phastCons_max,GERP_mean,GERP_max" > "$output_file"

    tail -n +2 "$initial_data"  | while IFS=',' read -r _ _ chr start end _; do

        #### PhyloP (pp) values: 241-way mammalian alignment ####
        echo "zoonomia_mean=$( $bigWigSummary_exe -type=mean $zoonomia_phylo_bw $chr $start $end 1 2>&1 )" >> errors.log
        echo "zoonomia_max=$(   $bigWigSummary_exe -type=max $zoonomia_phylo_bw $chr $start $end 1 2>&1 )" >> errors.log
        
        zoonomia_mean=$( $bigWigSummary_exe -type=mean $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
        zoonomia_max=$(   $bigWigSummary_exe -type=max $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
        
        
        ##### phastCons (pc) values: 100way-alignment #####
        echo "max_pc =   $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1" >> errors.log
        echo "mean_pc = $bigWigSummary_exe -type=mean $phast_bw $chr $start $end 1 2>&1" >> errors.log
        
        mean_pc=$( $bigWigSummary_exe -type=mean $phast_bw $chr $start $end 1 2>&1 )
        max_pc=$(   $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1 )
	    

        #### GERP Scores: 100way-alignment ####
        echo "gerp_mean = $bigWigSummary_exe -type=mean $gerp_bw $chr $start $end 1 2>&1" >> errors.log
        echo "gerp_max =   $bigWigSummary_exe -type=max $gerp_bw $chr $start $end 1 2>&1" >> errors.log

        gerp_mean=$( $bigWigSummary_exe -type=mean $gerp_bw $chr $start $end 1 2>&1 )
        gerp_max=$(   $bigWigSummary_exe -type=max $gerp_bw $chr $start $end 1 2>&1 )


        ## Convert any missing data to NA
        missing_value() {

            local var="$1"

            test_count=$(echo "$var" | wc -w)
            
            if [ "$test_count" -ne 1 ]; then var="NA"; fi
        
            echo "$var"
        }

        mean_pc=$(missing_value "$mean_pc")
        max_pc=$(missing_value "$max_pc")
        zoonomia_mean=$(missing_value "$zoonomia_mean")
        zoonomia_max=$(missing_value "$zoonomia_max")
        gerp_mean=$(missing_value "$gerp_mean")
        gerp_max=$(missing_value "$gerp_max")

        ## Values to output file 
        echo "$zoonomia_mean,$zoonomia_max,$mean_pc,$max_pc,$gerp_mean,$gerp_max" >> "$output_file"

    done

fi



