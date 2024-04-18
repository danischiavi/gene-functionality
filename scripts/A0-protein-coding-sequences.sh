#!/bin/bash
#
# Script Name: protein-coding-sequences.sh
#
# Author: Daniela Schiavinato
# Last edited: April 2024
#
# Description: 

#   - Generate Functional datasets: 
#                                 functional-protein-exon2 & functional-protein-exon3 
#
#   - Extracts coordinates to generate negative control dataset downstream using negative-control.sh:                     
#                                 protein-coords-negative-control
#                               
# 
###########################################################################################################################
#### GENERAL SETUP ####
###########################################################################################################################

##### Input files ##### 
genome_annotations=$1
genome_seq=$2
protein_coding_refseq=$3

mkdir -p data/datasets

#### Final Output files #####
# Functional datasets 
protein_exon_two='data/datasets/functional-protein-exon2-dataset.csv'
protein_exon_three='data/datasets/functional-protein-exon3-dataset.csv'
# Coordenates for negative-control generation 
coding_negative_control='data/datasets/protein-coords-negative-control.csv'

#### Temporary ####
protein_exon_two_unsorted='data/datasets/functional-protein-exon2-unsorted.csv'
protein_exon_three_unsorted='data/datasets/functional-protein-exon3-unsorted.csv'
coding_negative_control_unsorted='data/datasets/protein-coords-negative-unsorted.csv'

##### Constrains ##### 
sample_size=10
lower_limit_protein='61'
upper_limit_protein='272'

###########################################################################################################################
#### FUNCTION DECLARATION ####
###########################################################################################################################

###########################
# Get info from the Human genome database gff file to filter out pi-coding sequences 
gff2Info() {                                                                                

    local exons=$1
    local genome_seq=$2

    coords_one=$(   awk 'NR==1 {print $1, $4, $5}' "$exons")                                                      # Required to generate the upstream negative control sequences
    coords_two=$(   awk 'NR==2 {print $1, $4, $5}' "$exons")                                                      # Exon two coordinates
    coords_three=$( awk 'NR==3 {print $1, $4, $5}' "$exons")                                                      # Exon three coordinates
    final_end=$(tail -1 "$exons" | awk '{print $1, $4, $5}')                                                      # Required to generate the downstream negative control sequences

    chr=$( echo "$coords_one" | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 )                     # Chromosome variable
    test=$( echo "$chr" | cut -c1 )                                                                               # Records any zeros in the chromosome variable
    other=$( echo "$chr" | cut -c2 )                                                                              # If zero is in chromosome variable, only record the single digit (ie: 01 becomes 1).
    mt_test=$( echo "$coords_one" | cut -d ' ' -f 1 )                                                             # Variable to check if gene is located on the mitochondrial genome.

    # Reformat chr variable or rename to allow it to be filtered out
    if [ -z "$chr" ]                                                                                             # If chromosome variable empty (genes/mRNA that have been removed)
    then
        chr=26   
    elif [[ "$mt_test" == "NC_012920.1" ]]                                                                       # If gene is encoded on the mitochondrial genome
    then
        chr=25     
    elif [[ "$test" == "0" ]]                                                                                    # If chromosome variable begins with zero, then rename as a single digit (ie: 01 becomes 1).
    then
        chr="$other"    
    fi

    # Process exon data to create a dataset
    if [ "$chr" -le '23' ]; then  
    
        IFS=' ' read -r _ start_one _           <<< "$coords_one"                                                 # Starting and end coordinates of exons
        IFS=' ' read -r _ start_two end_two     <<< "$coords_two"
        IFS=' ' read -r _ start_three end_three <<< "$coords_three"
        
        if [[ "$chr" == 23 ]]; then chr=X; fi                                                                     # Chromosome X is NC_000023, but should be recorded as X in the final dataset for readability.
        
        seq_two=$(   grep -w "chromosome $chr" "$genome_seq" | cut -f 2 | cut -c"$start_two"-"$end_two" )         # Exons sequences 
        seq_three=$( grep -w "chromosome $chr" "$genome_seq" | cut -f 2 | cut -c"$start_three"-"$end_three" )       
        end_final=$( echo "$final_end" | cut -d ' ' -f 3 )                                                        # End position of final exon
	
        len_two=$(( end_two - start_two ))                                                                        # Length of exons
        len_three=$(( end_three - start_three ))          

        # Exclude empty and with any unknown nucleotides (N) sequences
        if [ -n "$seq_two" ] && [ -n "$seq_three" ] && [[ "$seq_two" != *"N"* ]] && [[ "$seq_three" != *"N"* ]]; then
       
	    # Exclude sequences out of length limits 
	        if { [ "$len_two" -ge "$lower_limit_protein" ] && [ "$len_two" -le "$upper_limit_protein" ]; } && \
            { [ "$len_three" -ge "$lower_limit_protein" ] && [ "$len_three" -le "$upper_limit_protein" ]; }; then 
                
                selected_ids+=("$random_id")                            
                protein_count=$(echo "${#selected_ids[@]}") 

                echo "RNA$protein_count,Yes,chr$chr,$start_two,$end_two,$seq_two" >> "$protein_exon_two_unsorted"
                echo "RNA$protein_count,Yes,chr$chr,$start_three,$end_three,$seq_three" >> "$protein_exon_three_unsorted"

                # Reverse transcripts can alter order of start/end positions
                if [ "$start_one" -gt "$end_final" ]; then                              # Negative strand                                                     
                    
                    distance_exons=$(( start_two - end_three ))
		            # To generate negative control sequences that are the same length as exons two and three
                    echo "chr$chr,$end_final,$start_one,$len_two,$len_three,$distance_exons" >> "$coding_negative_control_unsorted"
                
                else

                    distance_exons=$(( start_three - end_two ))        
                    echo "chr$chr,$start_one,$end_final,$len_two,$len_three,$distance_exons" >> "$coding_negative_control_unsorted"
                
                fi
            fi 
	    fi
	fi
}

###########################             

###########################################################################################################################
#### EXTRACT SEQUENCES ####
###########################################################################################################################

if [ ! -s "$protein_exon_two" ] || [ ! -s "$protein_exon_three" ]; then

    echo "ID,Functional,Chromosome,Start_gene,End_gene,Sequence" > "$protein_exon_two_unsorted"
    echo "ID,Functional,Chromosome,Start_gene,End_gene,Sequence" > "$protein_exon_three_unsorted"
    echo "Chromosome,Start,End,Length_exon2,Length_exon3,Distance_between_exons" > "$coding_negative_control_unsorted"

    declare -a "selected_ids=()"
    protein_count=0

    while [ "$protein_count" -lt "$sample_size" ]; do                                               # D: count in the gff2Info function --> it would be better to have it here --> TO THINK

        random_id=$(shuf -n 1 "$protein_coding_refseq")

        if [[ ! " ${selected_ids[@]} " =~ " $random_id " ]]; then

            grep "exon-$random_id" "$genome_annotations" > data/exons                               # Grep annotation from Reference Genome (NCBI) according to protein-coding genes (HGNC)
        
            if [ "$(wc -l < data/exons)" -ge 4 ]; then gff2Info data/exons "$genome_seq"; fi        # At least 4 exons 
            
        fi
    done 


    # Sort sequences and remove temporary files 
    sort_functional(){

        local unsorted_file=$1
        local output_file=$2

        awk -F ',' 'NR > 1 {print $1}' "$unsorted_file" > protein-id-column
        awk -F ',' 'NR > 1 {print $2 "," $3 "," $4 "," $5 "," $6}' "$unsorted_file" | sort -t ',' -k2,2 -k3,3n -k4,4n > protein-sorted_columns
        (echo "ID,Functional,Chromosome,Start,End,Sequence"; paste -d ',' protein-id-column protein-sorted_columns) > "$output_file"

        rm -rf protein-id-column
        rm -rf protein-sorted_columns
        rm -rf "$unsorted_file"

    }

    sort_functional "$protein_exon_two_unsorted" "$protein_exon_two"
    sort_functional "$protein_exon_three_unsorted" "$protein_exon_three"

    # Sort coordenates for negative control 
    awk 'NR == 1 {print $0; next} {print $0 | "sort -t, -k1,1 -k2,2n -k3,3n"}' "$coding_negative_control_unsorted" > "$coding_negative_control"

    rm -rf "$coding_negative_control_unsorted"

fi
