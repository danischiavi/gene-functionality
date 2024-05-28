#!/bin/bash
#
# Script Name: non-coding-sequences.sh
#
# Author: Daniela Schiavinato
# Last edited: April 2024
#
# Description: 

#   - Generate Functional datasets: 
#                                 functional-short-ncrna, 
#                                 functional-lncrna-exon1, functional-lncrna-exon2, 
#
#   - Extracts coordinates to generate negative control dataset downstream using negative-control.sh:                     
#                                 lncrna-coords-negative-control,
#                                 short-ncrna-coords-negative-control
# 

###########################################################################################################################
#### GENERAL SETUP ####
###########################################################################################################################

##### Input files ##### 
rnacentral_coords=$1
rnacentral_short_ncrna_seqs=$2
rnacentral_pre_mirna_seqs=$3
genome_seq=$4
interaction_database=$5

#### Final Output files #####
mkdir -p data/datasets

# Functional datasets 
short_ncrna='data/datasets/functional-short-ncrna-dataset.csv'
# Coordenates for negative-control generation 
short_negative_control='data/datasets/short-ncrna-coords-negative-control.csv'

##### Temporary files ####
int_database_file='data/datasets/interaction-database-tmp'
short_info='data/datasets/short-info'
short_tmp='data/datasets/short-tmp'
short_ncrna_unsorted='data/datasets/functional-short-ncrna-unsorted.csv'
short_negative_control_unsorted='data/datasets/short-ncrna-coords-negative-unsorted.csv'

##### Constrains ##### 
sample_size=100                     # Number of sequences for each type of RNA 
lower_limit_short='71'               # Given by the size distribution analysis (10&90% percentile)
upper_limit_short='142'

## Interaction database ## 
awk '/^>/ {print substr($1, 2)}' "$interaction_database" > "$int_database_file"
declare -a "IDs_interaction=()" 
mapfile -t IDs_interaction < "$int_database_file"

###########################################################################################################################
#### FUNCTION DECLARATION ####
###########################################################################################################################

# Set variables to filter out ncRNA sequences from the RNAcentral database
set_variables() {                                         
    
    local id=$1
    
    meta=$( grep -m 1 "$id" "$rnacentral_coords" )
    strand=$( echo "$meta" | cut -f 6 )

    IFS=$'\t ' read -r chr zero_start end _ <<< "$meta"                                                 # zero_start: 0-start bed format  

    # Remove if no sequence available     
    if [ "$chr" != 'chrM' ] && [ "$chr" != 'chrY' ]; then                                  
                                                
        len=$(( end-zero_start ))  
            
        if [ "$len" -ge "$lower_limit_short" ] && [ "$len" -le "$upper_limit_short" ]; then  
           
            # add 1 to change coordinate from 0-start
            start=$((zero_start + 1))  

            # Add random id to selected_ids arrray to keep track                                                                  
            selected_ids+=("$random_id")
            short_count="${#selected_ids[@]}"

            # Output functional dataset
            echo -e "$chr\t$zero_start\t$end\tRNA$short_count\t.\t$strand" >> "$short_info" 
                                                 

            if [ "$start" -gt "$end" ]; then                                                            # Reverse transcripts can alter order of start/end positions
                
                echo "$chr,$end,$start,$len,$strand" >> "$short_negative_control_unsorted"

            else

                echo "$chr,$start,$end,$len,$strand" >> "$short_negative_control_unsorted"

            fi
        fi 
    fi
}

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

# Populate arrays with non-coding-sequences IDs column for searching within all ncrna RNAcentral database 
declare -a "IDs_ncrna=()"  
mapfile -t IDs_ncrna < <(cut -f 4 -d $'\t' "$rnacentral_coords")


if [ ! -s "$short_ncrna" ]; then 

	if [ ! -s "$short_info" ]; then
    
		declare -a "IDs_short_ncrna=()"
    	mapfile -t IDs_short_ncrna < <(cut -f 1 -d $'\t' "$rnacentral_short_ncrna_seqs" | awk '{print $1}')
    
    	declare -a "IDs_pre_mirna=()"
    	mapfile -t IDs_pre_mirna < <(cut -f 1 -d $'\t' "$rnacentral_pre_mirna_seqs" | awk '{print $1}')

    	declare -a "selected_ids=()"
    	short_count=0

    	while [ "$short_count" -lt "$sample_size" ]; do

        	short_leftover="${#IDs_short_ncrna[@]}"
    
        	if [ "$short_leftover" -gt 0 ]; then

            	random_id=$(printf "%s\n" "${IDs_short_ncrna[@]}" | shuf -n 1)

            	# Find index of random id
            	for index in "${!IDs_short_ncrna[@]}"; do                                                     
                	if [[ "${IDs_short_ncrna[$index]}" == "$random_id" ]]; then
                    	del_index=$index
                    	break
                	fi
            	done

            	# Remove random id from short_ncrna array to avoid repeats and make downstream selection faster 
            	IDs_short_ncrna=( "${IDs_short_ncrna[@]:0:$del_index}" "${IDs_short_ncrna[@]:$((del_index + 1))}" )                                   

				if [[ ! " ${IDs_interaction[@]} " =~ " $random_id " ]]; then
            	
					if [[ "${IDs_ncrna[@]}" =~ "$random_id" ]]; then

                		set_variables "$random_id" 
                     
            		fi 
				fi

        	else 

            	random_id=$(printf "%s\n" "${IDs_pre_mirna[@]}" | shuf -n 1)

				if [[ ! " ${IDs_interaction[@]} " =~ " $random_id " ]]; then

            		if [[ ! " ${selected_ids[@]} " =~ " $random_id " ]]; then

                		if [[ "${IDs_ncrna[@]}" =~ "$random_id" ]]; then

                    		set_variables "$random_id" 
        
                		fi
					fi
            	fi
        	fi    
    	done
	fi
	
	## Extract sequences from genome ##
	if [ ! -s "$short_tmp" ]; then

		bedtools getfasta -fi "$genome_seq" -bed "$short_info" -fo "$short_tmp" -s -name -tab  
			# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
			# -name Use the name field and coordinates for the FASTA header
			# tab Report extract sequences in a tab-delimited format instead of in FASTA format.
	
		
	
	fi

	## Format Final file ## 
	if [ ! -s "$shor_ncrna_unsorted" ]; then
	
		reformat_file "$short_tmp" "$short_ncrna_unsorted"	
		sort_output "$short_ncrna_unsorted" "$short_ncrna" 

	fi
	
    # Sort coordenates for negative control 
	{
		echo "Chromosome,Start,End,Length,Strand"
    	sort -t, -k1,1 -k2,2n -k3,3n "$short_negative_control_unsorted"
	} > "$short_negative_control"

	rm -rf rm -rf "$short_info"
	rm -rf rm -rf "$short_tmp"
	rm -rf rm -rf "$short_ncrna_unsorted"
    rm -rf "$short_negative_control_unsorted" 

fi

rm -rf "$int_database_file"