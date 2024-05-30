#!/bin/bash
#
# Script Name: non-coding-sequences.sh
#
# Author: Daniela Schiavinato
# Last edited: April 2024
#
# Description: 

#   - Generate Functional datasets: 
#                                 functional-lncrna-exon1 & functional-lncrna-exon2
#
#   - Extracts coordinates to generate negative control dataset downstream using negative-control.sh:                     
#                                 lncrna-coords-negative-control

###########################################################################################################################
#### GENERAL SETUP ####
###########################################################################################################################

##### Input files ##### 
rnacentral_coords=$1
rnacentral_lncrna_seqs=$2
genome_seq=$3
interaction_database=$4

#### Final Output files #####
mkdir -p data/datasets

# Functional datasets 
lncrna_exon_one='data/datasets/functional-lncrna-exon1-dataset.csv' 
lncrna_exon_two='data/datasets/functional-lncrna-exon2-dataset.csv' 
# Coordenates for negative-control generation 
lncrna_negative_control='data/datasets/lncrna-coords-negative-control.csv'

#### Temporary files ####
int_database_file='data/datasets/interaction-database-tmp'
lncrna_exon_one_info='data/datasets/lncrna-exon-one-info'
lncrna_exon_two_info='data/datasets/lncrna-exon-two-info'
lncrna_exon_one_tmp='data/datasets/lncrna-exon-one-tmp'
lncrna_exon_two_tmp='data/datasets/lncrna-exon-two-tmp'

lncrna_exon_one_unsorted='data/datasets/functional-lncrna-exon1-unsorted.csv' 
lncrna_exon_two_unsorted='data/datasets/functional-lncrna-exon2-unsorted.csv'  
lncrna_negative_control_unsorted='data/datasets/lncrna-coords-negative-unsorted.csv'

##### Constrains ##### 
sample_size=1000                     # Number of sequences for each type of RNA 
lower_limit_lncrna='74'              # Given by the size distribution analysis (10&90% percentile)
upper_limit_lncrna='1500'

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

	if [ -n "$meta" ]; then

    	chr=$(  echo "$meta" | cut -f 1 )
		strand=$( echo "$meta" | cut -f 6 )
                             
    	if [ "$chr" != 'chrM' ] && [ "$chr" != 'chrY' ]; then 

        	exon_count=$( echo "$meta" | cut -f 10 )

        	if [ "$exon_count" -ge 2 ]; then                                                                # Filter for Multiexonic lncrna
                        
        		# Length of the required exons
            	len_first=$(  echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $1}')               # Length exons within range
            	len_last=$( echo "$meta" | awk -F'\t' '{print $11}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')
				
            	# Coordinates 
                seq_start=$(echo "$meta" | awk -F'\t' '{print $2}' )                                	# 0-start
                seq_end=$(echo "$meta" | awk -F'\t' '{print $3}' )    

				relative_start_first=$(  echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $1}') # Position relative to seq_start extracted from (bed format) for exons. For exon one is zero     
                start_first=$(( seq_start + relative_start_first ))                                  
            	end_first=$((              start_first + len_first ))                             

                relative_start_last=$( echo "$meta" | awk -F'\t' '{print $12}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')  # Extracted this way since just need the coordinates (not seq)
                start_last=$(( seq_start + relative_start_last )) 
                end_last=$((             start_last + len_last ))

				# On the positive strand, exons are named 1,2,3 from left to right, while the negative is right to left 
				# The order the len and position of the exons is described on the bed file is from left to right for both strands

					if [ "$strand" == '+' ]; then 
						
						len_one="$len_first"
						start_one="$start_first" 
						end_one="$end_first" 

						len_two=$(echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $2}') 
						relative_start_two=$( echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $2}')
                		start_two=$(( seq_start + relative_start_two ))                                       
                		end_two=$((              start_two + len_two ))
            
			    	elif [ "$strand" == '-' ]; then 

						len_one="$len_last"
						start_one="$start_last" 
						end_one="$end_last" 

						len_two=$( echo "$meta" | awk -F'\t' '{print $11}' | awk -v exon_count="$(( exon_count - 1 ))" -F',' '{print $exon_count}') 
						relative_start_two=$( echo "$meta" | awk -F'\t' '{print $12}' | awk -v exon_count="$(( exon_count - 1 ))" -F',' '{print $2}')
                		start_two=$(( seq_start + relative_start_two ))                                       
                		end_two=$((              start_two + len_two ))

					fi

					if { [ "$len_one" -ge "$lower_limit_lncrna" ] && [ "$len_one" -le "$upper_limit_lncrna" ]; } && { [ "$len_two" -ge "$lower_limit_lncrna" ] && [ "$len_two" -le "$upper_limit_lncrna" ]; }; then  
                
				    	(( lncrna_count ++ )) # RNA count

						echo -e "$chr\t$start_one\t$end_one\tRNA$lncrna_count.$id\t.\t$strand" >> "$lncrna_exon_one_info"
            			echo -e "$chr\t$start_two\t$end_two\tRNA$lncrna_count.$id\t.\t$strand" >> "$lncrna_exon_two_info" 
          	                                                       
                    	echo "$chr,$seq_start,$seq_end,$len_one,$len_two,$strand" >> "$lncrna_negative_control_unsorted"      # To generate negative control sequences that are the same length as exons two and three      	
            	fi
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
		split(parts[1], header_parts, ".")
        gene_id = header_parts[2]
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

find_index() {

    local value="$1"
    shift
    local IDs_lncrna=("$@")
    
    for i in "${!IDs_lncrna[@]}"; do
        if [[ "${IDs_lncrna[$i]}" == "$value" ]]; then
            echo "$i"
            return
        fi
    done
    
    # Return a non-zero status if the value is not found
    return 1
}
###########################################################################################################################
#### EXTRACT SEQUENCES ####
###########################################################################################################################

if [ ! -s "$lncrna_exon_one" ] || [ ! -s "$lncrna_exon_two" ]; then

	if [ ! -s "$lncrna_exon_one_info" ] || [ ! -s "$lncrna_exon_two_info" ]; then  
    
		declare -a "IDs_lncrna=()"
    	mapfile -t IDs_lncrna < <(cut -f 1 -d $'\t' "$rnacentral_lncrna_seqs" | awk '{print $1}')
    
    	declare -a "selected_ids=()"   
		declare -a "possible_ids=()"                                                                                 # Keeps track of selected random IDs
    	
		lncrna_count=0
		
		# Remove ids from curated interaction database to exclude them on the functional dataset
		for id in "${IDs_interaction[@]}"; do 

			index=$(find_index "$id" "${IDs_lncrna[@]}")
			unset 'IDs_lncrna[index]'
        	IDs_lncrna=("${IDs_lncrna[@]}")
		
		done

    	while [ "$lncrna_count" -lt "$sample_size" ]; do

			if [ "${#IDs_lncrna[@]}" -gt 0 ]; then
        
				random_id=$(printf "%s\n" "${IDs_lncrna[@]}" | shuf -n 1)

            	set_variables "$random_id" 

				# Remove random id from array not to be chosen again 
				index=$(find_index "$random_id" "${IDs_lncrna[@]}")
				unset 'IDs_lncrna[index]'
        		IDs_lncrna=("${IDs_lncrna[@]}")
			else
				echo "not enogh functional sequences for selected sample size"
				break  
			fi
    	done
	fi

	## Extract sequences from genome ##
	if [ ! -s "$lncrna_exon_one_tmp" ]; then  
		bedtools getfasta -fi "$genome_seq" -bed "$lncrna_exon_one_info" -fo "$lncrna_exon_one_tmp" -s -name -tab  
			# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
			# -name Use the name field and coordinates for the FASTA header
			# tab Report extract sequences in a tab-delimited format instead of in FASTA format.
	fi 

	if [ ! -s "$lncrna_exon_two_tmp" ]; then 
		bedtools getfasta -fi "$genome_seq" -bed "$lncrna_exon_two_info" -fo "$lncrna_exon_two_tmp" -s -name -tab 
	fi

	## Format Final file ## 
	if [ ! -s "$lncrna_exon_one_unsorted" ]; then 
		reformat_file "$lncrna_exon_one_tmp" "$lncrna_exon_one_unsorted"
		sort_output "$lncrna_exon_one_unsorted" "$lncrna_exon_one"
	fi

	if [ ! -s "$lncrna_exon_two_unsorted" ]; then 
		reformat_file "$lncrna_exon_two_tmp" "$lncrna_exon_two_unsorted"
		sort_output "$lncrna_exon_two_unsorted" "$lncrna_exon_two"
	fi

    
    # Sort coordenates for negative control 
	{
		echo "Chromosome,Start,End,Length_exon1,Length_exon2,Strand"
    	sort -t, -k1,1 -k2,2n -k3,3n "$lncrna_negative_control_unsorted" 
	} > "$lncrna_negative_control"

	rm -rf "$lncrna_exon_one_info" "$lncrna_exon_two_info"
	rm -rf "$lncrna_exon_one_tmp" "$lncrna_exon_two_tmp"
	rm -rf "$lncrna_exon_one_unsorted" "$lncrna_exon_two_unsorted"		
    rm -rf "$lncrna_negative_control_unsorted"

fi

rm -rf "$int_database_file"


