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
protein_exon_two_info='data/datasets/protein-exon-two-info'
protein_exon_three_info='data/datasets/protein-exon-three-info'
protein_exon_two_tmp='data/datasets/protein-exon-two-tmp'
protein_exon_three_tmp='data/datasets/protein-exon-three-tmp'
protein_exon_two_unsorted='data/datasets/functional-protein-exon2-unsorted.csv'
protein_exon_three_unsorted='data/datasets/functional-protein-exon3-unsorted.csv'
coding_negative_control_unsorted='data/datasets/protein-coords-negative-unsorted.csv'

##### Constrains ##### 
sample_size=1000
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

    coords_one=$(   awk 'NR==1 {print $1, $4, $5, $7}' "$exons")                                                      # Required to generate the upstream negative control sequences
    coords_two=$(   awk 'NR==2 {print $1, $4, $5}' "$exons")                                                      # Exon two coordinates
    coords_three=$( awk 'NR==3 {print $1, $4, $5}' "$exons")                                                      # Exon three coordinates
    final_end=$(tail -1 "$exons" | awk '{print $1, $4, $5}')                                                      # Required to generate the downstream negative control sequences
	
	strand=$( echo "$coords_one" | cut -d ' ' -f 4 )

    chr=$( echo "$coords_one" | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 )                     # Chromosome variable
    test=$( echo "$chr" | cut -c1 )                                                                               # Records any zeros in the chromosome variable
    other=$( echo "$chr" | cut -c2 )                                                                              # If zero is in chromosome variable, only record the single digit (ie: 01 becomes 1).
    mt_test=$( echo "$coords_one" | cut -d ' ' -f 1 )                                                             # Variable to check if gene is located on the mitochondrial genome.

    # Reformat chr variable or rename to allow it to be filtered out
    if [ -z "$chr" ]; then                                                                                        # If chromosome variable empty (genes/mRNA that have been removed)
        chr=26   
    elif [[ "$mt_test" == "NC_012920.1" ]]; then                                                                 # If gene is encoded on the mitochondrial genome
        chr=25     
    elif [[ "$test" == "0" ]]; then                                                                              # If chromosome variable begins with zero, then rename as a single digit (ie: 01 becomes 1).
        chr="$other"    
    fi

    # Process exon data to create a dataset
    if [ "$chr" -le '23' ]; then  
    
        IFS=' ' read -r _ start_one _           <<< "$coords_one"                                                 # Starting and end coordinates of exons
        IFS=' ' read -r _ start_two end_two     <<< "$coords_two"
        IFS=' ' read -r _ start_three end_three <<< "$coords_three"

	    if [[ "$chr" == 23 ]]; then chr=X; fi                                                                     # Chromosome X is NC_000023, but should be recorded as X in the final dataset for readability.
        
        end_final=$( echo "$final_end" | cut -d ' ' -f 3 )                                                        # End position of final exon
		start_two_zero=$(( start_two - 1 ))
		start_three_zero=$(( start_three -1 ))

        len_two=$(( end_two - start_two ))                                                                        # Length of exons
        len_three=$(( end_three - start_three ))          
       
	    # Exclude sequences out of length limits 
	    if { [ "$len_two" -ge "$lower_limit_protein" ] && [ "$len_two" -le "$upper_limit_protein" ]; } && \
        { [ "$len_three" -ge "$lower_limit_protein" ] && [ "$len_three" -le "$upper_limit_protein" ]; }; then 
                
            selected_ids+=("$random_id")                            
            protein_count=$(echo "${#selected_ids[@]}") 

            echo -e "chr$chr\t$start_two_zero\t$end_two\tRNA$protein_count\t.\t$strand" >> "$protein_exon_two_info"
            echo -e "chr$chr\t$start_three_zero\t$end_three\tRNA$protein_count\t.\t$strand" >> "$protein_exon_three_info"

            # Reverse transcripts can alter order of start/end positions
            if [ "$start_one" -gt "$end_final" ]; then                              # Negative strand                                                     
                    
		        # To generate negative control sequences that are the same length as exons two and three
                echo "chr$chr,$end_final,$start_one,$len_two,$len_three,$strand" >> "$coding_negative_control_unsorted"
                
            else
       
                echo "chr$chr,$start_one,$end_final,$len_two,$len_three,$strand" >> "$coding_negative_control_unsorted"
                
            fi
        fi 
	fi

}

###########################             
reformat_file() {

    local input_file=$1
    local output_file=$2

    awk -F'\t' '
    function calc_ambiguity(seq,   i, amb_nucleotides, total_nucleotides, ambiguity_percent) {
        amb_nucleotides = 0
        total_nucleotides = length(seq)
        
        for (i = 1; i <= total_nucleotides; i++) {
            if (index("RYMKSWBDHVNrymkswbdhvn", substr(seq, i, 1))) {
                amb_nucleotides++
            }
        }
        
        ambiguity_percent = (amb_nucleotides * 100) / total_nucleotides
        return ambiguity_percent
    	
		}

    {
        split($1, parts, "::")
        split(parts[2], coords, ":")
        chromosome = coords[1]
        split(coords[2], range, "-")
        start = range[1] + 1 
        end = substr(range[2], 1, index(range[2], "(") - 1)

        sequence = $2
        ambiguity = calc_ambiguity(sequence)

        if (ambiguity <= 5) {
            printf("Yes,%s,%s,%s,%s\n", chromosome, start, end, sequence)
        }
    }' "$input_file" >> "$output_file"

    rm -rf "$input_file"

}

###########################

sort_output(){

    local input_file=$1
    local output_file=$2

	local file_id=data/datasets/"$(basename "${input_file%.*}" | sed 's/-unsorted.csv//')" 

	sort -t ',' -k2,2 -k3,3n -k4,4n "$input_file" > "$file_id"-sorted-columns

	awk -v start=1 -v end="$sample_size" '
    BEGIN {
        for (i = start; i <= end; i++) {
            print "RNA" i
        }
    }' > "$file_id"-id-column

    (echo "ID,Functional,Chromosome,Start,End,Sequence"; paste -d ',' "$file_id"-id-column "$file_id"-sorted-columns) > "$output_file"

	rm -rf "$file_id"-id-column "$file_id"-sorted-columns
}

###########################################################################################################################
#### EXTRACT SEQUENCES ####
###########################################################################################################################

if [ ! -s "$protein_exon_two" ] || [ ! -s "$protein_exon_three" ]; then

	if [ ! -s "$protein_exon_two_info" ] || [ ! -s "$protein_exon_three_info" ]; then
    	
		declare -a "selected_ids=()"
    	protein_count=0

    	while [ "$protein_count" -lt "$sample_size" ]; do 
	                                              
        	random_id=$(shuf -n 1 "$protein_coding_refseq")

        	if [[ ! " ${selected_ids[@]} " =~ " $random_id " ]]; then

            	grep "exon-$random_id" "$genome_annotations" > data/exons                               # Grep annotation from Reference Genome (NCBI) according to protein-coding genes (HGNC)
        
            	if [ "$(wc -l < data/exons)" -ge 4 ]; then gff2Info data/exons "$genome_seq"; fi        # At least 4 exons 
            
        	fi
    	done 
	fi

	## Extract sequences from genome ##
	bedtools getfasta -fi "$genome_seq" -bed "$protein_exon_two_info" -fo "$protein_exon_two_tmp" -s -name -tab  
	bedtools getfasta -fi "$genome_seq" -bed "$protein_exon_three_info" -fo "$protein_exon_three_tmp" -s -name -tab 
	# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
	# -name Use the name field and coordinates for the FASTA header
	# tab Report extract sequences in a tab-delimited format instead of in FASTA format.

	reformat_file "$protein_exon_two_tmp" "$protein_exon_two_unsorted"
	reformat_file "$protein_exon_three_tmp" "$protein_exon_three_unsorted"

    sort_output "$protein_exon_two_unsorted" "$protein_exon_two"
    sort_output "$protein_exon_three_unsorted" "$protein_exon_three"


    # Sort coordenates for negative control 
	{
		echo "Chromosome,Start,End,Length_exon2,Length_exon3,Strand"
    	sort -t, -k1,1 -k2,2n -k3,3n "$coding_negative_control_unsorted" 
	} > "$coding_negative_control"

	## Remove excess files ## 
	rm -rf "$protein_exon_two_info" "$protein_exon_three_info"
	rm -rf "$protein_exon_two_tmp" "$protein_exon_three_tmp"
	rm -rf "$protein_exon_two_unsorted" "$protein_exon_three_unsorted"
    rm -rf "$coding_negative_control_unsorted"

fi
