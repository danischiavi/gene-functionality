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

#### Final Output files #####
# Functional datasets 
lncrna_exon_one='data/functional-lncrna-exon1-dataset.csv' 
lncrna_exon_two='data/functional-lncrna-exon2-dataset.csv' 
# Coordenates for negative-control generation 
lncrna_negative_control='data/lncrna-coords-negative-control.csv'

##### Constrains ##### 
sample_size=1000                     # Number of sequences for each type of RNA 
lower_limit_lncrna='74'              # Given by the size distribution analysis (10&90% percentile)
upper_limit_lncrna='1149'

###########################################################################################################################
#### FUNCTION DECLARATION ####
###########################################################################################################################

# Set variables to filter out ncRNA sequences from the RNAcentral database
set_variables() {                                         
    
    local id=$1
    
    seq=$(  grep -m 1 "$id" "$rnacentral_lncrna_seqs" | cut -f 2 )
    meta=$( grep -m 1 "$id" "$rnacentral_coords" )
    chr=$(  echo "$meta" | cut -f 1 )
                                         
    if [ -n "$seq" ] && [ "$chr" != 'chrM' ] && [ "$chr" != 'chrY' ]; then 

        exon_count=$( echo "$meta" | cut -f 10 )

        if [ "$exon_count" -ge 2 ]; then                                                                # Filter for Multiexonic lncrna
                        
        # Length of the required exons
            len_one=$(  echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $1}')               # Length exons within range
            len_two=$(echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $2}') 
            len_last=$( echo "$meta" | awk -F'\t' '{print $11}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')

            if { [ "$len_one" -ge "$lower_limit_lncrna" ] && [ "$len_one" -le "$upper_limit_lncrna" ]; } && { [ "$len_two" -ge "$lower_limit_lncrna" ] && [ "$len_two" -le "$upper_limit_lncrna" ]; }; then  
                
            # Coordinates to extract sequence from RNAcentral
                seq_start_zero=$(echo "$meta" | awk -F'\t' '{print $2}' )                                # 0-start. Function in progress for this (see ./bin/draft)
                seq_start=$((seq_start_zero + 1))

                relative_start_one=$(  echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $1}') # Position relative to seq_start extracted from (bed format) for exons. It's actually zero     
                start_index_one=$((relative_start_one + 1))                                              # index to extract sequece: +1 to account for the first relative start pos being 0
                end_index_one=$((relative_start_one + len_one))
                seq_one=$( echo "$seq" | cut -c $start_index_one-$end_index_one)                         # Sequence of exon 1 extracted from RNAcentral (fasta file converted to csv)

                start_index_two=$((end_index_one + 1))                                                   # Field 12 are the coordinates relative to seq_start but considering the introns,     # which are not present in the downloaded RNAcentral lncrna sequence -> to extract exon seq   # use values relative to exon 1                                          
                end_index_two=$((end_index_one + len_two))
                seq_two=$( echo "$seq" | cut -c $start_index_two-$end_index_two)
                        
            # Coordinates in genome 
                start_one=$((   $seq_start + $relative_start_one ))                                  
                end_one=$((                $start_one + $len_one - 1))                                   # -1 to make the end coordinate inclusive (to match UCSC browser)

                relative_start_two=$(  echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $2}')
                start_two=$((   $seq_start + $relative_start_two ))                                         
                end_two=$((                $start_two + $len_two - 1))

                relative_start_last=$( echo "$meta" | awk -F'\t' '{print $12}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')  # Extracted this way since just need the coordinates (not seq)
                start_last=$(( $seq_start + $relative_start_last )) 
                end_last=$((             $start_last + $len_last - 1))
            
            # RNA count
                selected_ids+=("$random_id")
                lncrna_count=$(echo "${#selected_ids[@]}")
            
                echo RNA$lncrna_count,Yes,"$chr","$start_one","$end_one","$seq_one" >> "$lncrna_exon_one"
                echo RNA$lncrna_count,Yes,"$chr","$start_two","$end_two","$seq_two" >> "$lncrna_exon_two"  

                if [ "$start_one" -gt "$end_last" ]; then                                                 # Reverse transcripts can alter order of start/end positions
                                                                             
                    echo "$chr,$end_last,$start_one,$len_one,$len_two" >> "$lncrna_negative_control"      # To generate negative control sequences that are the same length as exons two and three
                        
                else

                    echo "$chr,$start_one,$end_last,$len_one,$len_two" >> "$lncrna_negative_control"
                fi
            fi
        fi
    fi
}

###########################################################################################################################
#### EXTRACT SEQUENCES ####
###########################################################################################################################

# Populate arrays with non-coding-sequences IDs column for searching within all ncrna RNAcentral database 
declare -a "IDs_ncrna=()"  
mapfile -t IDs_ncrna < <(cut -f 4 -d $'\t' "$rnacentral_coords")

if [ ! -s "$lncrna_exon_one" ] || [ ! -s "$lncrna_exon_two" ]; then 

    echo "ID,Functional,Chromosome,Start,End,Sequence" >> "$lncrna_exon_one"
    echo "ID,Functional,Chromosome,Start,End,Sequence" >> "$lncrna_exon_two"
    echo "Chromosome,Start,End,Length_exon1,Length_exon2" >> "$lncrna_negative_control" 

    declare -a "IDs_lncrna=()"
    mapfile -t IDs_lncrna < <(cut -f 1 -d $'\t' "$rnacentral_lncrna_seqs" | awk '{print $1}')
    
    declare -a "selected_ids=()"                                                                                    # Keeps track of selected random IDs
    lncrna_count=0

    while [ "$lncrna_count" -lt "$sample_size" ]; do
    
        random_id=$(printf "%s\n" "${IDs_lncrna[@]}" | shuf -n 1)                                                   # Select a random ID from the lncrna list
                                             
        if [[ ! " ${selected_ids[@]} " =~ " $random_id " ]]; then                                                   # Select no repeated IDs
       
            if [[ "${IDs_ncrna[@]}" =~ "$random_id" ]]; then 
    
                set_variables "$random_id" 

            fi
        fi
    done
fi



