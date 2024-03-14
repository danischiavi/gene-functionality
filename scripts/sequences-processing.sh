#!/bin/bash
#
# Script Name: sequences-processing.sh
#
# Author: Daniela Schiavinato
# Last edited: March 2024
#
# Description: 

#   - Generate Functional datasets: 
#                                 functional-short-ncrna, 
#                                 functional-lncrna-exon1, functional-lncrna-exon2, 
#                                 functional-protein-exon2, functioanl-protein-exon3 
#
#   - Extracts coordinates to generate negative control dataset downstream using negative-control.sh:
#                                 protein-coords-negative-control
#                                 lncrna-coords-negative-control
#                                 short-ncrna-coords-negative-control
# 

###########################################################################################################################
#### GENERAL SETUP ####
###########################################################################################################################

# Input files # 
rnacentral_coords=$1
rnacentral_lncrna_seqs=$2
rnacentral_short_ncrna_seqs=$3
rnacentral_pre_mirna_seqs=$4
genome_annotations=$5
genome_seq=$6
protein_coding_refseq=$7

# Final Output files #
# Functional datasets 
lncrna_exon_one='data/functional-lncrna-exon1-dataset.csv' 
lncrna_exon_two='data/functional-lncrna-exon2-dataset.csv' 
short_ncrna='data/functional-short-ncrna-dataset.csv'
protein_exon_two='data/functional-protein-exon2-dataset.csv'
protein_exon_three='data/functional-protein-exon3-dataset.csv'

# Coordenates for negative-control generation 
coding_negative_control='data/protein-coords-negative-control.csv'
lncrna_negative_control='data/lncrna-coords-negative-control.csv'
short_negative_control='data/short-ncrna-coords-negative-control.csv'


# Constrains # 
sample_size=1000                     # Number of sequences for each type of RNA 
lower_limit_lncrna='74'              # Given by the size distribution analysis (10&90% percentile)
upper_limit_lncrna='1149'
lower_limit_short='71' 
upper_limit_short='142'
lower_limit_protein='61'
upper_limit_protein='272'


###########################################################################################################################
#### FUNCTION DECLARATION ####
###########################################################################################################################

###########################
# Set variables to filter out ncRNA sequences from the RNAcentral database
set_variables() {                                         
    
    local id=$1
    local seq_csv=$2
    
    seq=$(  grep -m 1 "$id" "$seq_csv" | cut -f 2 )
    meta=$( grep -m 1 "$id" "$rnacentral_coords" )
    chr=$(  echo "$meta" | cut -f 1 )
    status='pass'                                           # Set status to filter out sequences which are redifine as 'not-pass'
        
    if [ -z "$seq" ]; then                                  # Remove if no sequence available
        status='not-pass'
        return
	elif [ $chr == 'chrM' ] || [ $chr == 'chrY' ]; then     # Remove ncRNA from mitochondria and Y chromosome
        status='not-pass'
        return
    fi
}

###########################
# Get info from the Human genome database gff file to filter out pi-coding sequences 
gff2Info() {                                                                                

    local exons=$1
    local genome_seq=$2

    coords_one=$(   awk 'NR==1 {print $1, $4, $5}' "$exons")                                                    # Required to generate the upstream negative control sequences
    coords_two=$(   awk 'NR==2 {print $1, $4, $5}' "$exons")                                                    # Exon two coordinates
    coords_three=$( awk 'NR==3 {print $1, $4, $5}' "$exons")                                                    # Exon three coordinates
    final_end=$(tail -1 "$exons" | awk '{print $1, $4, $5}')                                                    # Required to generate the downstream negative control sequences

    chr=$( echo $coords_one | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 )                     # Chromosome variable
    test=$( echo $chr | cut -c1 )                                                                               # Records any zeros in the chromosome variable
    other=$( echo $chr | cut -c2 )                                                                              # If zero is in chromosome variable, only record the single digit (ie: 01 becomes 1).
    mt_test=$( echo $coords_one | cut -d ' ' -f 1 )                                                             # Variable to check if gene is located on the mitochondrial genome.

    # Reformat chr variable or rename to allow it to be filtered out
    if [ -z "$chr" ]                                                                                            # If chromosome variable empty (genes/mRNA that have been removed)
    then
        chr=26   
    elif [[ "$mt_test" == "NC_012920.1" ]]                                                                      # If gene is encoded on the mitochondrial genome
    then
        chr=25     
    elif [[ "$test" == "0" ]]                                                                                   # If chromosome variable begins with zero, then rename as a single digit (ie: 01 becomes 1).
    then
        chr="$other"    
    fi

    # Process exon data to create a dataset
    if [ "$chr" -le '23' ]; then  
    
        IFS=' ' read -r _ start_one end_one     <<< "$coords_one"                                               # Starting and end coordinates of exons
        IFS=' ' read -r _ start_two end_two     <<< "$coords_two"
        IFS=' ' read -r _ start_three end_three <<< "$coords_three"
        
        if [[ "$chr" == 23 ]]; then chr=X; fi                                                                   # Chromosome X is NC_000023, but should be recorded as X in the final dataset for readability.
        
        seq_two=$(   grep -w "chromosome $chr" "$genome_seq" | cut -f 2 | cut -c$start_two-$end_two )           # Exons sequences 
        seq_three=$( grep -w "chromosome $chr" "$genome_seq" | cut -f 2 | cut -c$start_three-$end_three )       
        end_final=$( echo $final_end | cut -d ' ' -f 3 )                                                        # End position of final exon
	
        len_two=$(( $end_two - $start_two ))                                                                    # Length of exons
        len_three=$(( $end_three - $start_three ))          
        
        # Exclude empty and with any unknown nucleotides (N) sequences
        if [ ! -z "$seq_two" ] && [ ! -z "$seq_three" ] && [[ "$seq_two" != *"N"* ]] && [[ "$seq_three" != *"N"* ]]; then
       
	    # Exclude sequences out of length limits 
	        if ([ "$len_two" -ge "$lower_limit_protein" ] && [ "$len_two" -le "$upper_limit_protein" ]) && ([ "$len_three" -ge "$lower_limit_protein" ] && [ "$len_three" -le "$upper_limit_protein" ]); then 
                
                selected_ids+=("$random_id")                            # D: probably is better if its in the while structure below instead
                protein_count=$(echo "${#selected_ids[@]}") 

                echo RNA$protein_count,Yes,chr$chr,$start_two,$end_two,$seq_two >> "$protein_exon_two"
                echo RNA$protein_count,Yes,chr$chr,$start_three,$end_three,$seq_three >> "$protein_exon_three"
    
                if [ "$start_one" -gt "$end_final" ]; then                                                       # Reverse transcripts can alter order of start/end positions
               
		            # To generate negative control sequences that are the same length as exons two and three
                    echo chr$chr,$end_final,$start_one,$len_two,$len_three >> "$coding_negative_control"
                
                else

                    echo chr$chr,$start_one,$end_final,$len_two,$len_three >> "$coding_negative_control"
                
                fi
            fi 
	    fi
	fi
}

###########################             


###########################################################################################################################
#### EXTRACT SEQUENCES ####
###########################################################################################################################

# Populate arrays with short/lncrna IDs column for searching within all ncrna RNAcentral database 
declare -a "IDs_ncrna=()"  
mapfile -t IDs_ncrna < <(cut -f 4 -d $'\t' "$rnacentral_coords")

############
## LNCRNA ##
############

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
    
                set_variables "$random_id" "$rnacentral_lncrna_seqs" 

                if [ "$status" != 'not-pass' ]; then                                                                # Status returned by function if empty sequence or mitocondrial/Ychr
         
                    exon_count=$( echo "$meta" | cut -f 10 )

                    if [ "$exon_count" -ge 2 ]; then                                                                # Filter for Multiexonic lncrna
                        
                        # Length of the required exons
                        len_one=$(  echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $1}')               # Length exons within range
                        len_two=$(echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $2}') 
                        len_last=$( echo "$meta" | awk -F'\t' '{print $11}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')

                        if ([ "$len_one" -ge "$lower_limit_lncrna" ] && [ "$len_one" -le "$upper_limit_lncrna" ]) && ([ "$len_two" -ge "$lower_limit_lncrna" ] && [ "$len_two" -le "$upper_limit_lncrna" ])  
                        then

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
            fi
        fi
    done
fi

#################
## SHORT-NCRNA ##
#################

if [ ! -s "$short_ncrna" ]; then 

    echo "ID,Functional,Chromosome,Start,End,Sequence" >> "$short_ncrna"
    echo "Chromosome,Start,End,Length" >> "$short_negative_control" 

    declare -a "IDs_short_ncrna=()"
    mapfile -t IDs_short_ncrna < <(cut -f 1 -d $'\t' "$rnacentral_short_ncrna_seqs" | awk '{print $1}')
    
    declare -a "IDs_pre_mirna=()"
    mapfile -t IDs_pre_mirna < <(cut -f 1 -d $'\t' "$rnacentral_pre_mirna_seqs" | awk '{print $1}')

    declare -a "selected_ids=()"
    short_count=0
    short_total="${#IDs_short_ncrna[@]}"

    while [ "$short_count" -lt "$sample_size" ]; do
    
        if [ "$short_count" -le "$short_total" ]; then

            random_id=$(printf "%s\n" "${IDs_short_ncrna[@]}" | shuf -n 1)                                   
        
            if [[ "${IDs_ncrna[@]}" =~ "$random_id" ]]; then

                set_variables "$random_id" "$rnacentral_short_ncrna_seqs" 

            fi    

        else 

            random_id=$(printf "%s\n" "${IDs_pre_mirna[@]}" | shuf -n 1)
                                           

            if [[ "${IDs_ncrna[@]}" =~ "$random_id" ]]; then

                set_variables "$random_id" "$rnacentral_pre_mirna_seqs"
        
            fi
        fi    

        if [ "$status" = 'pass' ]; then
        
            IFS=$'\t ' read -r chr zero_start end _ <<< "$meta"                                             # zero_start: 0-start bed format
            len=$(( $end-$zero_start ))  
            
            if [ "$len" -ge "$lower_limit_short" ] && [ "$len" -le "$upper_limit_short" ]; then  
           
                start=$((zero_start + 1))                                                                    # add 1 to change coordinate from 0-start
                selected_ids+=("$random_id")
                short_count="${#selected_ids[@]}"
                                                             
                echo RNA$short_count,Yes,"$chr","$start","$end","$seq" >> "$short_ncrna"
                
                if [ "$start" -gt "$end" ]; then                                                            # Reverse transcripts can alter order of start/end positions
                
                    echo "$chr,$end,$start,$len" >> "$short_negative_control"

                else

                    echo "$chr,$start,$end,$len" >> "$short_negative_control"

                fi
            fi 
        fi

    done
fi

########################
## PROTEIN-CODING-RNA ##
########################

if [ ! -s "$protein_exon_two" ] || [ ! -s "$protein_exon_three" ]; then

    echo "ID,Functional,Chromosome,Start,End,Sequence" >> "$protein_exon_two"
    echo "ID,Functional,Chromosome,Start,End,Sequence" >> "$protein_exon_three"
    echo "Chromosome,Start,End,Length_exon2,Length_exon3" >> "$coding_negative_control"

    declare -a "selected_ids=()"
    protein_count=0

    while [ "$protein_count" -lt "$sample_size" ]; do                                               # D: count in the gff2Info function --> it would be better to have it here --> TO THINK

        random_id=$(shuf -n 1 "$protein_coding_refseq")

        if [[ ! " ${selected_ids[@]} " =~ " $random_id " ]]; then

            grep "exon-$random_id" "$genome_annotations" > data/exons                               # Grep annotation from Reference Genome (NCBI) according to protein-coding genes (HGNC)
        
            if [ "$(wc -l < data/exons)" -ge 4 ]; then gff2Info data/exons "$genome_seq"; fi        # At least 4 exons 
            
        fi
    done 
fi
