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
sample_size=100                     # Number of sequences for each type of RNA 
lower_limit_lncrna='74'              # Given by the size distribution analysis (10&90% percentile)
upper_limit_lncrna='1149'

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
    chr=$(  echo "$meta" | cut -f 1 )
	strand=$( echo "$meta" | cut -f 6 )

                                         
    if [ "$chr" != 'chrM' ] && [ "$chr" != 'chrY' ]; then 

        exon_count=$( echo "$meta" | cut -f 10 )

        if [ "$exon_count" -ge 2 ]; then                                                                # Filter for Multiexonic lncrna
                        
        # Length of the required exons
            len_one=$(  echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $1}')               # Length exons within range
            len_two=$(echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $2}') 
            len_last=$( echo "$meta" | awk -F'\t' '{print $11}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')

            if { [ "$len_one" -ge "$lower_limit_lncrna" ] && [ "$len_one" -le "$upper_limit_lncrna" ]; } && { [ "$len_two" -ge "$lower_limit_lncrna" ] && [ "$len_two" -le "$upper_limit_lncrna" ]; }; then  
                
            # Coordinates to extract sequence from RNAcentral
                seq_start_zero=$(echo "$meta" | awk -F'\t' '{print $2}' )                                # 0-start
                seq_start=$((seq_start_zero + 1))
                relative_start_one=$(  echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $1}') # Position relative to seq_start extracted from (bed format) for exons. For exon one is zero     
                
            # Coordinates in genome 
                start_one=$((   seq_start + relative_start_one ))                                  
                start_one_zero=$(( start_one - 1 ))
				end_one=$((             start_one + len_one - 1))                                   	 # -1 to make the end coordinate inclusive (to match UCSC browser)

                relative_start_two=$(  echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $2}')
                start_two=$((   seq_start + relative_start_two ))
				start_two_zero=$(( start_two - 1 ))                                         
                end_two=$((             start_two + len_two - 1))

                relative_start_last=$( echo "$meta" | awk -F'\t' '{print $12}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')  # Extracted this way since just need the coordinates (not seq)
                start_last=$(( seq_start + relative_start_last )) 
                end_last=$((          start_last + len_last - 1))
            
            # RNA count
                selected_ids+=("$random_id")
                lncrna_count=$(echo "${#selected_ids[@]}")
	
				echo -e "$chr\t$start_one_zero\t$end_one\tRNA$lncrna_count\t.\t$strand" >> "$lncrna_exon_one_info"
            	echo -e "$chr\t$start_two_zero\t$end_two\tRNA$lncrna_count\t.\t$strand" >> "$lncrna_exon_two_info" 

                if [ "$start_one" -gt "$end_last" ]; then                                                			# Reverse transcripts can alter order of start/end positions
                                                                             
                    echo "$chr,$end_last,$start_one,$len_one,$len_two,$strand" >> "$lncrna_negative_control_unsorted"      # To generate negative control sequences that are the same length as exons two and three
                        
                else

                    echo "$chr,$start_one,$end_last,$len_one,$len_two,$strand" >> "$lncrna_negative_control_unsorted"
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
        chromosome = coords[1]
        split(coords[2], range, "-")
        start = range[1] + 1 
        end = substr(range[2], 1, index(range[2], "(") - 1)

        sequence = $2
        ambiguity = calc_ambiguity(sequence)

        if (ambiguity <= 5) {
            printf("%s,%s,%s,%s\n", chromosome, start, end, sequence)
        }
    }' "$input_file" >> "$output_file"

    rm -rf "$input_file"

}

###########################

sort_output(){

    local file_with_duplicates=$1
    local output_file=$2

    local file_id="$(basename "${file_with_duplicates%.*}" | sed 's/-dataset//')" 
    local file_no_duplicates_unsorted="$file_id"-unsorted

	# Remove duplicates 
    awk -F',' 'BEGIN { OFS="\t" }
            {
                key = $1 FS $2 FS $3 FS $4
                count[key]++
                if (count[key] <= 1)
                    print $0
            }' "$file_with_duplicates" > "$file_no_duplicates_unsorted"


	sort -t ',' -k1,1 -k2,2n -k3,3n "$file_no_duplicates_unsorted" > "$file_id"-sorted-columns

    numb_seqs=$(wc -l < "$file_id"-sorted-columns)

	awk -v start="$initial_data_seq_numb" -v end="$(( initial_data_seq_numb + numb_seqs - 1 ))" '
    BEGIN {
        for (i = start; i <= end; i++) {
            print "RNA" i
        }
    }' > "$file_id"-id-column

    (echo "ID,Functional,Chromosome,Start,End,Sequence"; paste -d ',' "$file_id"-id-column "$file_id"-sorted-columns) > "$output_file"

}

###########################################################################################################################
#### EXTRACT SEQUENCES ####
###########################################################################################################################

# Populate arrays with non-coding-sequences IDs column for searching within all ncrna RNAcentral database 
declare -a "IDs_ncrna=()"  
mapfile -t IDs_ncrna < <(cut -f 4 -d $'\t' "$rnacentral_coords")

if [ ! -s "$lncrna_exon_one" ] || [ ! -s "$lncrna_exon_two" ]; then 

    echo "ID,Functional,Chromosome,Start,End,Sequence" >> "$lncrna_exon_one_unsorted"
    echo "ID,Functional,Chromosome,Start,End,Sequence" >> "$lncrna_exon_two_unsorted"
    echo "Chromosome,Start,End,Length_exon1,Length_exon2" >> "$lncrna_negative_control_unsorted" 

    declare -a "IDs_lncrna=()"
    mapfile -t IDs_lncrna < <(cut -f 1 -d $'\t' "$rnacentral_lncrna_seqs" | awk '{print $1}')
    
    declare -a "selected_ids=()"                                                                                    # Keeps track of selected random IDs
    lncrna_count=0

    while [ "$lncrna_count" -lt "$sample_size" ]; do
    
        random_id=$(printf "%s\n" "${IDs_lncrna[@]}" | shuf -n 1)                                                   # Select a random ID from the lncrna list

		if [[ ! " ${IDs_interaction[@]} " =~ " $random_id " ]]; then

        	if [[ ! " ${selected_ids[@]} " =~ " $random_id " ]]; then                                                   # Select no repeated IDs
       
            	if [[ "${IDs_ncrna[@]}" =~ "$random_id" ]]; then 
    
                	set_variables "$random_id" 

            	fi
        	fi
		fi
		
    done

	## Extract sequences from genome ##
	bedtools getfasta -fi "$genome_seq" -bed "$lncrna_exon_one_info" -fo "$lncrna_exon_one_tmp" -s -name -tab  
	bedtools getfasta -fi "$genome_seq" -bed "$lncrna_exon_two_info" -fo "$lncrna_exon_two_tmp" -s -name -tab 
	# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
	# -name Use the name field and coordinates for the FASTA header
	# tab Report extract sequences in a tab-delimited format instead of in FASTA format.
	#rm -rf "$lncrna_exon_one_info"
	#rm -rf rm -rf "$lncrna_exon_two_info"
	
	## Format Final file ## 
	reformat_file "$lncrna_exon_one_tmp" "$lncrna_exon_one_unsorted"
	reformat_file "$lncrna_exon_two_tmp" "$lncrna_exon_two_unsorted"

    sort_functional "$lncrna_exon_one_unsorted" "$lncrna_exon_one"
    sort_functional "$lncrna_exon_two_unsorted" "$lncrna_exon_two"

    # Sort coordenates for negative control 
    awk 'NR == 1 {print $0; next} {print $0 | "sort -t, -k1,1 -k2,2n -k3,3n"}' "$lncrna_negative_control_unsorted" > "$lncrna_negative_control"

    rm -rf "$lncrna_negative_control_unsorted"

fi



