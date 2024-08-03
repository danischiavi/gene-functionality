#!/bin/bash
#
# Script Name: 
#
# Author: Daniela Schiavinato
# Last edited: 
#

#                               
# 
###########################################################################################################################
#### GENERAL SETUP ####
###########################################################################################################################

##### Input files ##### 
genome_annotations=$1
genome_seq=$2
protein_coding_refseq=$3

dir=data/datasets/testModel
mkdir -p "$dir"

#### Final Output files #####
# Functional datasets 
protein_exon_two="${dir}/testModel-protein-exon2-dataset.csv"
protein_exon_three="${dir}/testModel-protein-exon3-dataset.csv"

#### Temporary ####
protein_exon_two_info="${dir}/testModel-protein-exon2-info"
protein_exon_three_info="${dir}/testModel-protein-exon3-info"
protein_exon_two_tmp="${dir}/testModel-protein-exon2-tmp"
protein_exon_three_tmp="${dir}/testModel-protein-exon3-tmp"
protein_exon_two_unsorted="${dir}/testModel-protein-exon2-unsorted"
protein_exon_three_unsorted="${dir}/testModel-protein-exon3-unsorted"

##### Constrains ##### 
sample_size=500
lower_limit_protein='61'
upper_limit_protein='272'

###########################################################################################################################
#### FUNCTION DECLARATION ####
###########################################################################################################################

###########################
# Get info from the Human genome database gff file to filter out pi-coding sequences 
gff2Info() {                                                                                

    local exons=$1
   
    coords_two=$(   awk 'NR==2 {print $1, $4, $5, $7}' "$exons")                                                      # Exon two coordinates
    coords_three=$( awk 'NR==3 {print $1, $4, $5}' "$exons")                                                      # Exon three coordinates
    
	strand=$( echo "$coords_two" | cut -d ' ' -f 4 )

    chr=$( echo "$coords_two" | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 )                     # Chromosome variable
    test=$( echo "$chr" | cut -c1 )                                                                               # Records any zeros in the chromosome variable
    other=$( echo "$chr" | cut -c2 )                                                                              # If zero is in chromosome variable, only record the single digit (ie: 01 becomes 1).
    mt_test=$( echo "$coords_two" | cut -d ' ' -f 1 )                                                             # Variable to check if gene is located on the mitochondrial genome.

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
                                         
        IFS=' ' read -r _ start_two end_two _   <<< "$coords_two"
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

            echo -e "chr$chr\t$start_two_zero\t$end_two\tRNA$protein_count.$random_id\t.\t$strand" >> "$protein_exon_two_info"
            echo -e "chr$chr\t$start_three_zero\t$end_three\tRNA$protein_count.$random_id\t.\t$strand" >> "$protein_exon_three_info"

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
		split(parts[1], header_parts, ".")
        gene_id = header_parts[2]
        split(parts[2], coords, ":")
        chromosome = coords[1]
        split(coords[2], range, "-")
        start = range[1] + 1 
        end = substr(range[2], 1, index(range[2], "(") - 1)

        sequence = $2
        ambiguity = calc_ambiguity(sequence)

        if (ambiguity <= 5) {
            printf("Yes,%s,%s,%s,%s,%s\n", chromosome, start, end, sequence, gene_id)
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

    (echo "ID,Functional,Chromosome,Start,End,Sequence,GeneID"; paste -d ',' "$file_id"-id-column "$file_id"-sorted-columns) > "$output_file"

	rm -rf "$file_id"-id-column "$file_id"-sorted-columns
}

###########################################################################################################################
#### EXTRACT SEQUENCES ####
###########################################################################################################################

if [ ! -s "$protein_exon_two" ] || [ ! -s "$protein_exon_three" ]; then

	if [ ! -s "$protein_exon_two_info" ] || [ ! -s "$protein_exon_three_info" ]; then
    	
		# Get GeneIDs that has already been used to train/test predicting model 
		declare -a "used_ids=()"
		values=$(awk -F, '$3 == "Yes" {print $1}' results/protein-exon2-features-matrix)
		IFS=$'\n' read -r -d '' -a used_ids <<< "$values"

		# To avoid repeats
		declare -a "selected_ids=()"
    	
    	protein_count=0

    	while [ "$protein_count" -lt "$sample_size" ]; do 
	                                              
        	random_id=$(shuf -n 1 "$protein_coding_refseq")

			if [[ ! " ${selected_ids[@]} " =~ " $random_id " ]]; then

        		if [[ ! " ${used_ids[@]} " =~ " $random_id " ]]; then

            		grep "exon-$random_id" "$genome_annotations" > data/exons                               # Grep annotation from Reference Genome (NCBI) according to protein-coding genes (HGNC)
        
            		if [ "$(wc -l < data/exons)" -ge 4 ]; then gff2Info data/exons "$genome_seq"; fi        # At least 4 exons 
            
				fi
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

	## Remove excess files ## 
	rm -rf "$protein_exon_two_info" "$protein_exon_three_info"
	rm -rf "$protein_exon_two_tmp" "$protein_exon_three_tmp"
	rm -rf "$protein_exon_two_unsorted" "$protein_exon_three_unsorted"

fi


