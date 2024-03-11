#!/bin/bash
#
# Script Name: seq-conservation-features.sh
#
# Description: This script calculates the PhyloP, PhyloCons and GERP for for protein-coding exons, lncRNA exons, short ncRNAs 
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
phylo_bw=$3
phast_bw=$4
zoonomia_phylo_bw=$5
gerp_bw=$6

output_directory=data/conservation
mkdir -p "$output_directory"

output_file="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')conservation.csv"                  


#### Variables and databases #### --> D: still not sure which is the best way to handle this inputs . Some are relative paths to the directory of the project and other are absolute 



#### Extract conservation features for each set of chromosome coordinates #### 
if [ ! -s "$output_file" ]; then

    echo "241w_PP_mean,241w_PP_max,100w_PP_mean,100w_PP_max,100w_PC_mean,100w_PC_max,GERP_mean,GERP_max" > "$output_file"

    tail -n +2 "$initial_data"  | while IFS=',' read -r _ _ chr start end _; do

        ##### phastCons (pc) and phyloP (pp) values: 100way-alignment #####
	    echo "mean_pp = $bigWigSummary_exe -type=mean $phylo_bw $chr $start $end 1 2>&1" >> errors.log
        mean_pp=$( $bigWigSummary_exe -type=mean $phylo_bw $chr $start $end 1 2>&1 )
        max_pp=$(   $bigWigSummary_exe -type=max $phylo_bw $chr $start $end 1 2>&1 )
        mean_pc=$( $bigWigSummary_exe -type=mean $phast_bw $chr $start $end 1 2>&1 )
        max_pc=$(   $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1 )
	    echo "max_pc = $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1" >> errors.log
	
        #### PhyloP (pp) values: 241-way mammalian alignment ####
        zoonomia_mean=$( $bigWigSummary_exe -type=mean $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
        zoonomia_max=$(   $bigWigSummary_exe -type=max $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
    
        #### Obtain GERP Scores: 100w-alignment ####
        gerp_mean=$( $bigWigSummary_exe -type=mean $gerp_bw $chr $start $end 1 2>&1 )
        gerp_max=$(   $bigWigSummary_exe -type=max $gerp_bw $chr $start $end 1 2>&1 )


        ## Convert any missing data to NA
        missing_value() {

            local var="$1"

            test_count=$(echo "$var" | wc -w)
            
            if [ "$test_count" -ne 1 ]; then var="NA"; fi
        
            echo "$var"
        }

        mean_pp=$(missing_value "$mean_pp")
        max_pp=$(missing_value "$max_pp")
        mean_pc=$(missing_value "$mean_pc")
        max_pc=$(missing_value "$max_pc")
        zoonomia_mean=$(missing_value "$zoonomia_mean")
        zoonomia_max=$(missing_value "$zoonomia_max")
        gerp_mean=$(missing_value "$gerp_mean")
        gerp_max=$(missing_value "$gerp_max")

        ## Values to output file 
        echo "$zoonomia_mean,$zoonomia_max,$mean_pp,$max_pp,$mean_pc,$max_pc,$gerp_mean,$gerp_max" >> "$output_file"

    done

fi



