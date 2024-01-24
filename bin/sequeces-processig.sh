#!/bin/bash
#
# Script Name: slicing-database.sh
#
# Author: Daniela Schiavinato
# Last edited: November 2023
#
# Description: This script filters the functional long and short-ncRNAs from RNAcentral database and protein-coding-RNA from xxx
#


###########################################################################################################################

# General setup

###########################################################################################################################
pre_mirna_fasta=
short_ncrna_fasta=
lncrna_fasta=


### Input Files

# check if files exist 
if [ ! -f "$lncrna_fasta" ] || [ ! -f "$short_ncrna_fasta" ] || [ ! -f "$chr_coord_ncrna" ]  ## ADD THE OTHER FILES!
then
    echo "Initial data files do not exist"
    break 
else
    :
fi

## Convert input from FASTA to CSV file for easier parsability
fasta_formatter -i "$lncrna_fasta" -o "data/raw/rnacentral-lncrna-seq.csv" -t              # downloaded fasta file from RNAcentral with: lncRNA sequences
fasta_formatter -i "$short_ncrna_fasta" -o "data/raw/rnacentral-short-ncrna-seq.csv" -t    #                                             short-ncRNA sequences
fasta_formatter -i "$pre_mirna_fasta" -o "data/raw/rnacentral-pre-mirna-seq.csv" -t        #                                             precursor miRNA sequences

## Variables
lncrna_fasta='data/raw/rnacentral-lncrna.fasta'
short_ncrna_fasta='data/raw/rnacentral-short-ncrna.fasta'
pre_mirna_fasta='data/raw/rnacentral-pre-miRNA.fasta'

rnacentral_lncrna_seqs='data/raw/rnacentral-lncrna-seq.csv'
rnacentral_short_ncrna_seqs='data/raw/rnacentral-short-ncrna-seq.csv'
rnacentral_pre_mirna_seqs='data/raw/rnacentral-pre-mirna-seq.csv'

rnacentral_coords='data/raw/rnacentral_GRCh38_coordinates.bed'
RefSeq_protein_coding='data/raw/hgnc-protein-coding-RefSeq.txt'
genome_annotations='data/raw/GRCh38_p14_genomic.gff'
genome_csv='data/raw/GRCh38.p14_genome.csv'

### Output files

# File creation
file_creation() {

    local name=$1
    echo "ID,Functional,Chromosome,Start,End,Sequence" > data/$name-dataset.csv
}

names=("protein-exon2" "protein-exon3" "functional-lncrna-exon1" "functional-lncrna-exon2" "functional-short-ncrna")
for name in "${names[@]}"; do file_creation "$name"; done


# Set variables 
lncrna_exon_one='data/functional-lncrna-exon1-dataset.csv' 
lncrna_exon_two='data/functional-lncrna-exon2-dataset.csv' 
short_ncrna='data/functional-short-ncrna-dataset.csv'
protein_exon_two='data/protein-exon2-dataset.csv'
protein_exon_three='data/protein-exon3-dataset.csv'
negative_control='data/coords-for-negative-control.csv'


### Sequence length limits
lower_limit='75'
upper_limit='3000'
lower_limit_short='20' # CHECK

###########################################################################################################################

# NON-CODING RNA - Chromosome coordinates and sequence 

###########################################################################################################################

#### Function declaration

populate_array() {                                ## Populate arrays with short and lncrna ids for searching within all ncrna RNAcentral database
    
    local file=$1
    local name=$2
    local -n array_ref="IDs_$name"

    while IFS=$'\t' read -r ID_full _; do         # ID with description (rnacentral) 
       IFS=' ' read -r ID _     <<< "$ID_full"
       array_ref+=("$ID")
    done < "$file"
}


set_variables() {                                  ## Set variables to filter out ncRNA sequences from the RNAcentral database
    
    local id=$1
    local seq_csv=$2
    
    seq=$(  grep -m 1 "$id" "$seq_csv" | cut -f 2 )
    meta=$( grep -m 1 "$id" "$rnacentral_coords" )
    chr=$(  echo "$meta" | cut -f 1 )
    status='pass'                                  # Set status to filter out sequences which are redifine as 'not-pass'
        
    if [ -z "$seq" ]; then                               # Remove if no sequence available
        status='not-pass'
        return
	elif [ $chr == 'chrM' ] || [ $chr == 'chrY' ]; then  # Remove ncRNA from mitochondria and Y chromosome
        status='not-pass'
        return
    fi
}
                

## Arrays for searching short and lncrna ids within all ncrna RNAcentral database  

declare -a "IDs_ncrna=()"  
declare -a "IDs_lncrna=()"
declare -a "IDs_short_ncrna=()"
declare -a "IDs_pre_mirna=()"

# Populate lncRNA and short-ncRNA arrays
populate_array "$rnacentral_lncrna_seqs" "lncrna"
populate_array "$rnacentral_short_ncrna_seqs" "short_ncrna"
populate_array "$rnacentral_pre_mirna_seqs" "pre_mirna"

# Populate ncRNA array
while IFS=$'\t' read -r _ _ _ ID _; do 
    if [[ "$ID" != "${IDs_ncrna[@]}" ]]; then
        IDs_ncrna+=("$ID")
    else  
        :
    fi
done < "$rnacentral_coords"
   
## D: ADD CHECKS FOR ARRAYS!    


## LNCRNA - Extracting 1000 random sequences from RNAcentral database

lncrna_count=0
declare -a "selected_ids=()"                            # Keeps track of selected random IDs

while [ "$lncrna_count" -le 1000 ]; do
    
    random_id="${IDs_lncrna[RANDOM % ${#IDs_lncrna[@]}]}" # Select a random ID from the lncrna list
                                                     
    if [[ ! " ${selected_ids[@]} " =~ " $random_id " ]]; then  # Select no repeated IDs
       
        if [[ "${IDs_ncrna[@]}" =~ "$random_id" ]]; then 
    
            set_variables "$random_id" "$rnacentral_lncrna_seqs" 

            if [ "$status" != 'not-pass' ]; then                                                                      # Status returned by function if empty sequence or mitocondrial/Ychr
         
                exon_count=$( echo "$meta" | cut -f 10 )

                if [ "$exon_count" -ge 2 ]; then                                                                # Filter for Multiexonic lncrna
             
                    len_one=$(  echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $1}')               # Length exons within range
                    len_two=$(echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $2}') 
                    len_last=$( echo "$meta" | awk -F'\t' '{print $11}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')

                    if ([ "$len_one" -ge "$lower_limit" ] && [ "$len_one" -le "$upper_limit" ]) && ([ "$len_two" -ge "$lower_limit" ] && [ "$len_two" -le "$upper_limit" ])  
                    then
                        
                        seq_start=$(echo "$meta" | awk -F'\t' '{print $2}' )                                     # 0-start. Function in progress for this (see ./bin/draft)
                        relative_start_one=$(  echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $1}') # Position relative to Start (bed format) for exons     
                        relative_start_two=$(  echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $2}')     
                        relative_start_last=$( echo "$meta" | awk -F'\t' '{print $12}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')  
                    
                        start_one=$((   $seq_start + $relative_start_one + 1 ))                                  # +1 to account for the first relative start pos being 0
                        start_two=$((   $seq_start + $relative_start_two + 1 ))                                         
                        start_last=$(( $seq_start + $relative_start_last + 1 )) 
                        end_two=$((                    $start_two + $len_two ))
                        end_last=$((                 $start_last + $len_last ))
                        
                        selected_ids+=("$random_id")
                        lncrna_count=$(echo "${#selected_ids[@]}")
            
                        echo RNA$lncrna_count,Yes,"$chr","$start_one","$end_one","$seq"     >> "$lncrna_exon_one"
                        echo RNA$lncrna_count,Yes,"$chr","$start_two","$end_two","$seq" >> "$lncrna_exon_two"  

                        if [ "$start_one" -gt "$end_last" ]; then                                                       # Reverse transcripts can alter order of start/end positions
                                                                             
                            echo $chr,$end_last,$start_one,$len_one,$len_two >> "$negative_control"                   # To generate negative control sequences that are the same length as exons two and three
                        
                        else

                            echo $chr,$start_one,$end_last,$len_one,$len_two >> "$negative_control"
                        fi

                        
                    fi
                fi
            fi
        fi
    fi
done

## SHORT-NCRNA - Extracting 1000 random sequences from RNAcentral database

declare -a "selected_ids=()"
short_count=0
short_total=$(echo "${#IDs_short_ncrna[@]}")


while [ "$short_count" -le 1000 ]; do
    
    if [[ "$short_count" -le "short_total"]]; then

        random_id="${IDs_short_ncrna[RANDOM % ${#IDs_short_ncrna[@]}]}"                                                  # Select a random ID from the lncrna list

        if [[ "${IDs_ncrna[@]}" =~ "$random_id" ]]; then

            set_variables "$random_id" "$rnacentral_short_ncrna_seqs" 

        fi    

    else 

        random_id="${IDs_pre_mirna[RANDOM % ${#IDs_short_ncrna[@]}]}"                                                  # Select a random ID from the lncrna list

        if [[ "${IDs_ncrna[@]}" =~ "$random_id" ]]; then

            set_variables "$random_id" "$rnacentral_pre_mirna_seqs"
        
        fi
    fi    

    if [ "$status" = 'pass' ]; then
        
        IFS=$'\t ' read -r chr zero_start end _ <<< "$meta"                                             # zero_start: 0-start bed format
        len=$(( $end-$zero_start ))  
            
        if [ "$len" -ge "$lower_limit_short" ] && [ "$len" -le "$upper_limit" ]; then  
           
            start=$((zero_start+=1))                                                                    # add 1 to change coordinate from 0-start
            selected_ids+=("$random_id")
            short_count=$(echo "${#selected_ids[@]}") 
                                                             
            echo RNA$short_count,Yes,"$chr","$start","$end","$seq" >> "$short_ncrna"
                
            if [ "$start" -gt "$end" ]; then                                                            # Reverse transcripts can alter order of start/end positions
                
                echo $chr,$end,$start,$len,'NA' >> "$negative_control"

            else

                echo $chr,$start,$end,$len,'NA' >> "$negative_control"

            fi
        fi 
    fi
done


###########################################################################################################################

# PROTEIN CODING RNA - Chromosome coordinates and sequence 

###########################################################################################################################

##### Function declaration ####

gff2Info() {                                                                                ## Get info from the Human genome database gff file to filter out pi-coding sequences 

    local exons=$1
    local genome_csv=$2

    coords_one=$(   awk 'NR==1 {print $1, $4, $5}' "$exons")                                # Required to generate the upstream negative control sequences
    coords_two=$(   awk 'NR==2 {print $1, $4, $5}' "$exons")                                # Exon two coordinates
    coords_three=$( awk 'NR==3 {print $1, $4, $5}' "$exons")                                # Exon three coordinates
    final_end=$(tail -1 "$exons" | awk '{print $1, $4, $5}')                                # Required to generate the downstream negative control sequences

    chr=$( echo $coords_one | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 ) # Chromosome variable
    test=$( echo $chr | cut -c1 )                                                           # Records any zeros in the chromosome variable
    other=$( echo $chr | cut -c2 )                                                          # If zero is in chromosome variable, only record the single digit (ie: 01 becomes 1).
    mt_test=$( echo $coords_one | cut -d ' ' -f 1 )                                         # Variable to check if gene is located on the mitochondrial genome.

    # Reformat chr variable or rename to allow it to be filtered out
    if [ -z "$chr" ]                                                                        # If chromosome variable empty (genes/mRNA that have been removed)
    then
        chr=26   
    elif [[ "$mt_test" == "NC_012920.1" ]]                                                  # If gene is encoded on the mitochondrial genome
    then
        chr=25     
    elif [[ "$test" == "0" ]]                                                               # If chromosome variable begins with zero, then rename as a single digit (ie: 01 becomes 1).
    then
        chr="$other"    
    else
        :
    fi

    # Process exon data to create a dataset
    if [ "$chr" -le '23' ]  
    then
        IFS=' ' read -r _ start_one end_one     <<< "$coords_one"                                               # Starting and end coordinates of exons
        IFS=' ' read -r _ start_two end_two     <<< "$coords_two"
        IFS=' ' read -r _ start_three end_three <<< "$coords_three"
        
        if [[ "$chr" == 23 ]]; then chr=X; else :; fi                                                           # Chromosome X is NC_000023, but should be recorded as X in the final dataset for readability.
        
        seq_two=$(   grep -w "chromosome $chr" "$genome_csv" | cut -f 2 | cut -c$start_two-$end_two )           # Exons sequences 
        seq_three=$( grep -w "chromosome $chr" "$genome_csv" | cut -f 2 | cut -c$start_three-$end_three )       
        end_final=$( echo $final_end | cut -d ' ' -f 3 )                                                        # End position of final exon
	
        len_two=$(( $end_two - $start_two ))                                                                    # Length of exons
        len_three=$(( $end_three - $start_three ))          
        
        # Exclude empty and with any unknown nucleotides (N) sequences
        if [ ! -z "$seq_two" ] && [ ! -z "$seq_three" ] && [[ "$seq_two" != *"N"* ]] && [[ "$seq_three" != *"N"* ]]
        then
       
	    # Exclude sequences out of length limits 
	        if ([ "$len_two" -ge "$lower_limit" ] && [ "$len_two" -le "$upper_limit" ]) && ([ "$len_three" -ge "$lower_limit" ] && [ "$len_three" -le "$upper_limit" ]) 
	        then
                
                selected_ids+=("$random_id")                            # D: probably is better if its in the while structure below nstead
                protein_count=$(echo "${#selected_ids[@]}") 

                echo RNA$protein_count,Yes,chr$chr,$start_two,$end_two,$seq_two >> "$protein_exon_two"
                echo RNA$protein_count,Yes,chr$chr,$start_three,$end_three,$seq_three >> "$protein_exon_three"
    
                if [ "$start_one" -gt "$end_final" ]                                                            # Reverse transcripts can alter order of start/end positions
                then
		            # To generate negative control sequences that are the same length as exons two and three
                    echo chr$chr,$end_final,$start_one,$len_two,$len_three >> "$negative_control"
                else
                    echo chr$chr,$start_one,$end_final,$len_two,$len_three >> "$negative_control"
                fi
            fi 
	    fi
	fi
}

###############################

# Apply to each of 1000 random sequences

declare -a "selected_ids=()"
protein_count=0

while [ "$protein_count" -lt 1000 ]; do                             # D: count in the gff2Info function --> it would be better to have it here --> TO THINK

    random_id=$(shuf -n 1 "$RefSeq_protein_coding")

    if [[ ! " ${selected_ids[@]} " =~ " $random_id " ]]; then

        grep "exon-$random_id" "$genome_annotations" > data/exons        # Grep annotation from Reference Genome (NCBI) according to protein-coding genes (HGNC)
        
        if [ "$(wc -l < data/exons)" -ge 4 ]; then                # At least 4 exons 
            gff2Info data/exons "$genome_csv"
        fi
    fi
done 



###########################################################################################################################

# Generate functional FASTA file from CSV file

###########################################################################################################################

csv_to_fasta() {                     
    
    {
        while IFS=',' read -r id _ _ _ _ seq
        do
            echo -e ">$id\n$seq" >> data/$name-seq.fa
        
        done
    } < data/$name-dataset.csv   
}   

for name in "${names[@]}"; do csv_to_fasta "$name"; done
