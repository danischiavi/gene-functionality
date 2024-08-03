#!/bin/bash
#
# Script Name: 
#
# Author: Daniela Schiavinato
# Last edited: 
#
# Description: 



###########################################################################################################################
#### GENERAL SETUP ####
###########################################################################################################################

##### Input files ##### 
rnacentral_coords=$1
rnacentral_lncrna_seqs=$2
genome_seq=$3

#### Final Output files #####
dir="data/datasets/testModel"
mkdir -p "$dir"

# Functional datasets 
exon_one="${dir}/testModel-lncrna-exon1-dataset.csv"
exon_two="${dir}/testModel-lncrna-exon2-dataset.csv" 

#### Temporary files ####
exon_one_info="${dir}/testModel-lncrna-exon1-info"
exon_two_info="${dir}/testModel-lncrna-exon2-info"
exon_one_tmp="${dir}/testModel-lncrna-exon1-tmp"
exon_two_tmp="${dir}/testModel-lncrna-exon2-tmp"

exon_one_unsorted="${dir}/testModel-lncrna-exon1-unsorted.csv" 
exon_two_unsorted="${dir}/testModel-lncrna-exon2-unsorted.csv"  

##### Constrains ##### 
lower_limit='50'              
upper_limit='3000'
sample_size=500

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

				if { [ "$len_one" -ge "$lower_limit" ] && [ "$len_one" -le "$upper_limit" ]; } && { [ "$len_two" -ge "$lower_limit" ] && [ "$len_two" -le "$upper_limit" ]; }; then  
                
				    (( lncrna_count ++ )) # RNA count

					echo -e "$chr\t$start_one\t$end_one\tRNA$lncrna_count.$id\t.\t$strand" >> "$exon_one_info"
            		echo -e "$chr\t$start_two\t$end_two\tRNA$lncrna_count.$id\t.\t$strand" >> "$exon_two_info" 
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
        end = substr(range[2], 1, index(range[2], "(") -1)

        sequence = $2
        ambiguity = calc_ambiguity(sequence)

        if (ambiguity <= 5) {
            printf("Yes,%s,%s,%s,%s,%s\n", chromosome, start, end, sequence, gene_id)
        }
    }' "$input_file" >> "$output_file"

	# start + 1 since is following bed format 

    rm -rf "$input_file"
}

###########################

sort_output(){

    local input_file=$1
    local output_file=$2

	local file_id=data/datasets/"$(basename "${input_file%.*}" | sed 's/-unsorted.csv//')" 

	sort -t ',' -k2,2 -k3,3n -k4,4n "$input_file" > "$file_id"-sorted-columns

	awk -v start=1 -v end=$(( $(wc -l < "$input_file")))  '
    BEGIN {
        for (i = start; i <= end; i++) {
            print "RNA" i
        }
    }' > "$file_id"-id-column

    (echo "ID,Functional,Chromosome,Start,End,Sequence,GeneID"; paste -d ',' "$file_id"-id-column "$file_id"-sorted-columns) > "$output_file"

	rm -rf "$file_id"-id-column "$file_id"-sorted-columns
	rm -rf "$input_file"
}

###########################

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

if [ ! -s "$exon_one" ] || [ ! -s "$exon_two" ]; then

	if [ ! -s "$exon_one_info" ] || [ ! -s "$exon_two_info" ]; then  

		declare -a "IDs_lncrna=()"
    	mapfile -t IDs_lncrna < <(cut -f 1 -d $'\t' "$rnacentral_lncrna_seqs" | awk '{print $1}')
    
		# Get GeneIDs that has already been used to train/test predicting model 
		declare -a "used_ids=()"
		values=$(awk -F, '$3 == "Yes" {print $1}' results/lncrna-exon2-features-matrix)
		IFS=$'\n' read -r -d '' -a used_ids <<< "$values"

		#declare -a "testModel_IDs=()"
		#for id in "${IDs_lncrna[@]}"; do
		#	if [[ ! " ${used_ids[@]} " =~ " $id " ]]; then  

		#		testModel_IDs+=("$id")
		#	fi
		# done

		lncrna_count=0 

		while [ "$lncrna_count" -lt "$sample_size" ]; do

			if [ "${#IDs_lncrna[@]}" -gt 0 ]; then

				random_id=$(printf "%s\n" "${IDs_lncrna[@]}" | shuf -n 1)

				# Remove random id from array not to be chosen again
				index=$(find_index "$random_id" "${IDs_lncrna[@]}")
				unset 'IDs_lncrna[index]'
        		IDs_lncrna=("${IDs_lncrna[@]}")

				if [[ ! " ${used_ids[@]} " =~ " $random_id " ]]; then        
	
            	set_variables "$random_id"
			
				fi

			else
				echo "not enogh functional sequences for selected sample size"
				break  
			fi
		done 
	fi

	
	## Extract sequences from genome ##
	if [ ! -s "$exon_one_tmp" ]; then  
		bedtools getfasta -fi "$genome_seq" -bed "$exon_one_info" -fo "$exon_one_tmp" -s -name -tab  
			# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
			# -name Use the name field and coordinates for the FASTA header
			# tab Report extract sequences in a tab-delimited format instead of in FASTA format.
	fi 

	if [ ! -s "$lncrna_exon_two_tmp" ]; then 
		bedtools getfasta -fi "$genome_seq" -bed "$exon_two_info" -fo "$exon_two_tmp" -s -name -tab 
	fi

	## Format Final file ## 
	if [ ! -s "$exon_one_unsorted" ]; then 
		reformat_file "$exon_one_tmp" "$exon_one_unsorted"
		sort_output "$exon_one_unsorted" "$exon_one"
	fi

	if [ ! -s "$exon_two_unsorted" ]; then 
		reformat_file "$exon_two_tmp" "$exon_two_unsorted"
		sort_output "$exon_two_unsorted" "$exon_two"
	fi

else 
	echo "Files already exist" 
fi 
    

