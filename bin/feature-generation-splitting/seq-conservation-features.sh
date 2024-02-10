#!/bin/bash
#
# Script Name: seq-conservation-features.sh
#
# Description: This script calculates the PhyloP, PhyloCons and GERP for for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences.
#
# Input: $1 is the dataset 
#        $2 is the location of the folder with the local databases/datasets or version specific executables
# 
# Any additional required files, directories or dependencies will be requested when the script is run and require manual
# input.
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.
#
############################################################################################################################

# Calculating sequence conservation features

############################################################################################################################

initial_data=$1

output_directory=data/conservation
mkdir -p "$output_directory"

output_file="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')conservation.csv"                  # Define name and directory for output file
   

# Variables and databases
bigWigSummary_exe=bin/bigWigSummary
phylo_bw=/Volumes/archive/userdata/student_users/danielaschiavinato/dani-scratch/features-of-function-data/hg38.phyloP100way.bw
phast_bw=/Volumes/archive/userdata/student_users/danielaschiavinato/dani-scratch/features-of-function-data/hg38.phastCons100way.bw
zoonomia_phylo_bw=data/raw/hg38.phyloP241mammalian2020v2.bigWig
gerp_bw=/Volumes/archive/userdata/student_users/danielaschiavinato/dani-scratch/features-of-functional-human-genes/data/raw/gerp_conservation_scores.homo_sapiens.GRCh38.bw

if [ ! -e "$output_file" ]; then

    echo "241w_PP_mean,241w_PP_max,100w_PP_mean,100w_PP_max,100w_PC_mean,100w_PC_max,GERP_mean,GERP_max" > $output_file

    {
        read 
        while IFS=',' read -r _ _ chr start end _; do

            ######## Obtain phastCons (pc) and phyloP (pp) values for each set of chromosome coordinates
	        echo "mean_pp = $bigWigSummary_exe -type=mean $phylo_bw $chr $start $end 1 2>&1" >> errors.log
            mean_pp=$( $bigWigSummary_exe -type=mean $phylo_bw $chr $start $end 1 2>&1 )
            max_pp=$(   $bigWigSummary_exe -type=max $phylo_bw $chr $start $end 1 2>&1 )
            mean_pc=$( $bigWigSummary_exe -type=mean $phast_bw $chr $start $end 1 2>&1 )
            max_pc=$(   $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1 )
	        echo "max_pc = $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1" >> errors.log
	
            ######## Extract PhyloP (pp) values from 241-way mammalian alignment for each set of chromosome coordinates
            zoonomia_mean=$( $bigWigSummary_exe -type=mean $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
            zoonomia_max=$(   $bigWigSummary_exe -type=max $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
    
            ######## Obtain GERP Scores for each set of chromosome coordinates
            gerp_mean=$( $bigWigSummary_exe -type=mean $gerp_bw $chr $start $end 1 2>&1 )
            gerp_max=$(   $bigWigSummary_exe -type=max $gerp_bw $chr $start $end 1 2>&1 )

            ######## Convert any missing data to NAs (mean = mn, max = mx, phyloP = p, phastCons = c) --> D: declare a short function for this

            test_mnp=$( echo $mean_pp | wc -w )
            [[ "$test_mnp" -eq "1" ]] || mean_pp="NA"           # check and change the others if working 
	
            test_mxp=$( echo $max_pp | wc -w )
            if [[ "$test_mxp" -eq "1" ]]; then : ; else max_pp=NA ; fi
	
            test_mnc=$( echo $mean_pc | wc -w )
            if [[ "$test_mnc" -eq "1" ]]; then : ; else mean_pc=NA ; fi
	
            test_mxc=$( echo $max_pc | wc -w )
            if [[ "$test_mxc" -eq "1" ]]; then : ; else max_pc=NA; fi
	
            test_mnc=$( echo $zoonomia_mean | wc -w )
            if [[ "$test_mnc" -eq "1" ]]; then : ; else mean_pc=NA ; fi
	
            test_mxc=$( echo $zoonomia_max | wc -w )
            if [[ "$test_mxc" -eq "1" ]]; then : ; else max_pc=NA ; fi

            test_mnc=$( echo $gerp_mean | wc -w )
            if [[ "$test_mnc" -eq "1" ]]; then : ; else gerp_mean=NA ; fi
	
            test_mxc=$( echo $gerp_max | wc -w )
            if [[ "$test_mxc" -eq "1" ]]; then : ; else gerp_max=NA ; fi
    
        echo "$zoonomia_mean,$zoonomia_max,$mean_pp,$max_pp,$mean_pc,$max_pc,$gerp_mean,$gerp_max" >> "$output_file"

        done

    } < "$initial_data"

fi


