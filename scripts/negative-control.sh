#!/bin/bash
#
# Notes: 
#        initial_data: 
#        first refers to the first exon of the initial data, which not necessary corresponds to the exon number 1 
#        last refers to the second exon of the initial data, which not necessary corresponds to the exon number 2
#        single refers to sequences which have a single exon 
#
#
#### Variables and files ####
initial_data=$1
genome_seq=$2
genes_complement=$3

file_name=data/$(basename "${initial_data%.*}" | sed 's/-coords-negative-control//')                 

# Temporary files
first_negative_coords="$file_name"-first-negative-coords.csv
last_negative_coords="$file_name"-last-negative-coords.csv
single_negative_coords="$file_name"-negative-coords.csv

# Output files
single_negative_control="$file_name"-negative-control-dataset.csv   # File for short-ncrna is the final negative control dataset
first_negative_control="$file_name"-first-negative-control.csv      # File for lncrna and protein will be rename after this script to associated with the corresponding exon  
last_negative_control="$file_name"-last-negative-control.csv

###########################################################################################################################
#### FUNCTION DECLARATION ####
###########################################################################################################################

###########################
# To determinate percentage of ambiguous nucleotides on sequences (used to filter out if >5%)
ambiguity_percent() {
    
    seq=$1

    amb_nucleotides=$( echo "$seq" | grep -o '[RYMKSWBDHVNrymkswbdhvn]' | wc -l )
    total_nucleotides=$( echo -n "$seq" | wc -c )
    ambiguity_percent=$(echo "scale=0; ($amb_nucleotides * 100) / $total_nucleotides" | bc)

    echo "$ambiguity_percent"

}

###########################
# To extract coordinates and sequences of the negative control given a distance from the sequence

#### For sequences with more than 1 exon #### 
multiple_exons_negative_control() {

    local distance_to_seq=$1
                                                                                             
    distance_between_exons=$(( "$end_seq" - "$start_seq" - "$first_len" - "$last_len" ))                                                    

    # UPSTREAM
    last_exon_left_end=$((     "$start_seq" - "$distance_to_seq" ))
    last_exon_left_start=$(( "$last_exon_left_end" - "$last_len" ))                                                                         # length based on exons's length
    
    first_exon_left_end=$(( "$last_exon_left_start" - "$distance_between_exons" ))                              
    first_exon_left_start=$(( "$first_exon_left_end" - "$first_len" ))


    if [[ "$last_exon_left_end" -gt 0 ]] && [[ "$last_exon_left_start" -gt 0 ]]; then                                                                 # Filter out if negative coord
                             
        echo "$chr,$last_exon_left_start,$last_exon_left_end,$distance_to_seq,$last_len" >> "$last_negative_coords"

    fi

    if [[ "$first_exon_left_end" -gt 0 ]] && [[ "$first_exon_left_start" -gt 0 ]]; then 
            
        echo "$chr,$first_exon_left_start,$first_exon_left_end,$distance_to_seq,$first_len" >> "$first_negative_coords"

    fi
   

    ## DOWNSTREAM
    first_exon_right_start=$(( "$end_seq" + "$distance_to_seq" ))
    first_exon_right_end=$(( "$first_exon_right_start" + "$first_len" ))                                

    last_exon_right_start=$(( "$first_exon_right_end" + "$distance_between_exons" ))
    last_exon_right_end=$((  "$last_exon_right_start" + "$last_len" ))


    if [[ "$first_exon_right_start" -gt 0 ]] && [[ "$first_exon_right_end" -gt 0 ]]; then 
  
        echo "$chr,$first_exon_right_start,$first_exon_right_end,$distance_to_seq,$first_len" >> "$first_negative_coords"
    
    fi
    

    if [[ "$last_exon_right_start" -gt 0 ]] && [[ "$last_exon_right_end" -gt 0 ]]; then 
  
        echo "$chr,$last_exon_right_start,$last_exon_right_end,$distance_to_seq,$last_len" >> "$last_negative_coords"
        
    fi

}

#### For sequences with single exon (short-ncRNA) 
single_exon_negative_control() {

    local distance_to_seq=$1
    
    chromo=$( echo "$chr" | tr -d "chr" )

    # UPSTREAM
    left_end=$(( "$start_seq" - "$distance_to_seq" ))
    left_start=$(( "$left_end" - "$len" ))                                                              
    
    if [[ "$left_end" -gt 0 ]] && [[ "$left_start" -gt 0 ]]; then 
    
        echo "$chr,$left_start,$left_end,$distance_to_seq,$len" >> "$single_negative_coords"
    
    fi


    # DOWNSTREAM                       
    right_start=$(( "$end_seq" + "$distance_to_seq" ))
    right_end=$(( "$right_start" + "$len" ))                                
                              

    if [[ "$right_start" -gt 0 ]] && [[ "$right_end" -gt 0 ]]; then 

        echo "$chr,$right_start,$right_end,$distance_to_seq" >> "$single_negative_coords"
    
    fi

}

###########################
# To filter negative control sequences using RNAcentral and GENCODE

filter_out_functional(){
    
    local negative_coords=$1
    local exon=$2
    local negative_control=$3

    ## Temporary files
    closest_coords="$file_name"-"$exon"-closest.bed
    bed_coords_sorted="$file_name"-"$exon"-coords-sorted.bed

    # Format and sort coordinates for bedtools 
    awk -F',' 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' "$negative_coords" | sort -k1,1V -k2,2n | tr ' ' '\t' > "$bed_coords_sorted"
   
    # Find the closest gene complement region to negative control coordinates 
    bedtools closest -a "$bed_coords_sorted" -b "$genes_complement" -D ref | awk -F'\t' '!seen[$7,$8]++' > "$closest_coords"


    ## Format final output file
    while IFS=$'\t' read -r chr _ _ distance len _ start end closestdistance; do

        chromo=$( echo "$chr" | tr -d "chr" )

        if [ "$closestdistance" -ge 0 ]; then                   # if closest upstream initial negative control coord
        
            end=$(( "$start" + "$len" ))
        
        else                                                    # if closest downstream initial negative control coord

            start=$(( "$end" - "$len" ))
        
        fi
        
        # Calculate distance to gene border 
        distance_gene=$(( "$distance" + "$closestdistance" ))

        # Extract sequence from genome 
        seq=$( grep -w "chromosome $chromo" "$genome_seq" | cut -f 2 | cut -c"$start"-"$end" )

        # Remove negative control sequences that with an ambiguity percent > 5
        if [ -n "$seq" ] && [[ $( ambiguity_percent "$seq" ) -lt 5 ]]; then  
  
            (( count++ ))
            echo "RNA$count,No,$chr,$start,$end,$seq,$distance_gene" >> "$negative_control"
        
        fi
   
    done < "$closest_coords"

    #rm -rf "$bed_coords_sorted"
    
}

###########################

###########################################################################################################################
#### EXTRACT NEGATIVE CONTROLS ####
###########################################################################################################################

distances_to_seq=("1000" "10000" "100000" "1000000" "5000000")                             # Distances choosen for negative control (see manuscript for justification)

num_fields=$( awk -F',' '{print NF}' "$initial_data" | sort -nu | tail -n 1 )              # To define if it's single or multiple exon sequences


# Single exon sequences # 
if [[ "$num_fields" -eq 4 ]]; then                                                                      

    if [ ! -s "$single_negative_coords" ]; then 

        tail -n +2 "$initial_data"  | while IFS=, read -r chr start_seq end_seq len; do   

            for position in "${distances_to_seq[@]}"; do single_exon_negative_control "$position"; done

        done 

    fi

    
    if [ ! -s "$single_negative_control" ]; then

        echo "ID,Functional,Chromosome,Start,End,Sequence,Distance" > "$single_negative_control"
        count=$(( $(wc -l < "$initial_data") - 1 ))                                     # Start counting from last functional seq  
        filter_out_functional "$single_negative_coords" 'single' "$single_negative_control"

    fi

fi


# Multiple exons sequences #
if [[ "$num_fields" -eq 5 ]]; then 

    if [ ! -s "$first_negative_coords" ] || [ ! -s "$last_negative_coords" ]; then 

        tail -n +2 "$initial_data"  | while IFS=, read -r chr start_seq end_seq first_len last_len; do  

            for position in "${distances_to_seq[@]}"; do multiple_exons_negative_control "$position"; done

        done 

    fi


    if [ ! -s "$first_negative_control" ]; then

        echo "ID,Functional,Chromosome,Start,End,Sequence,Distance" > "$first_negative_control"
        count=$(( $(wc -l < "$initial_data") - 1 ))                                    
        filter_out_functional "$first_negative_coords" 'first' "$first_negative_control"

    fi


    if [ ! -s "$last_negative_control" ]; then

        echo "ID,Functional,Chromosome,Start,End,Sequence,Distance" > "$last_negative_control"
        count=$(( $(wc -l < "$initial_data") - 1 ))                                     
        filter_out_functional "$last_negative_coords" 'last' "$last_negative_control"

    fi

fi

