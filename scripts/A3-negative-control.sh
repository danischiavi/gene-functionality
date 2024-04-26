#!/bin/bash
# 
## This script generates the corresponding negative control dataset for the input file (protein-codingRNA, lncRNA or short-ncRNA)       
# - Outputs one file for each exon (Exon 2&3 for protein-coding; 1&2 for lncRNA and single exon for short-ncRNA)      
# - Note that A refers to the exon of a functional seq of the RNA type (Exon 1 for lncRNA and Exon 2 for protein-codingRNA);
# B to the second exon (Exon 2 for lncRNA and Exon 3 for protein-codingRNA); and
# single refers to sequences which have a single exon (short-ncRNA)
#  
######## General Idea ########
# 1) Reading line by line, define a sampling region from the gene to a fix value (first downstream) 
#    This fix value corresponds to the furthest away we want to sample a negative control (highest d - see line 15). 
# 2) Find the intersection between the sampling region and the "genes complement"
# 3) For each exon, find the closest region of the "genes complement" to a position=(Exon + d )  
# 4) Take a random sample within this closest region of length=(exon's length)
# 6) Repeat (3) & (4) for each d. Note values for d are pre-defined - see manuscript for justification 
# 7) Repeat (1)-(6) for upstream 
##############################

#### General Set Up ####

## Input files ##
# Of each genes from were functional seq was selected: chr, start coords of first exon, last coord of last exon, length exonA and length exonB  
initial_data=$1             
# csv file with genome sequences for each chromosome 
genome_seq=$2
# regions of genome lacking of genes (annotated from RNAcentral and Gencode)
genes_complement=$3

# distances from genome to find a region to sample the negative control. See manuscript for justification
distances_to_seq=("1000" "10000" "100000" "1000000" "5000000")                             

# Name for files
file_name=data/datasets/$(basename "${initial_data%.*}" | sed 's/-coords-negative-control//')                 

## Temporary files ##
first_negative_coords="$file_name"-exonA-negative-coords.csv
last_negative_coords="$file_name"-exonB-negative-coords.csv
single_negative_coords="$file_name"-negative-coords.csv

sampling_region="$file_name"-sampling-region
complement_sampling_region="$file_name"-complement-sampling-region
                                                                                                                                   
closest_input="$file_name"-closest-input
closest_output="$file_name"-closest-output
random_output="$file_name"-random-coords

## Final output files ##
# Note the output files depends on the initial_data input file 
negative_control_single="$file_name"-negative-control-dataset.csv   # File for short-ncrna is the final negative control dataset
negative_control_A="$file_name"-exonA-negative-control.csv          # File for lncrna and protein will be rename after this script to associated with the corresponding exon  
negative_control_B="$file_name"-exonB-negative-control.csv

###########################################################################################################################
#### FUNCTION DECLARATION ####
###########################################################################################################################
# Finds the intersection between the sampling region and the "genes complement", removing the regions which length is smaller than desired 
# $genes_complement is a global variable of the script which is assigned with a script input

intersect(){

    sampling_region=$1
    len=$2
    complement_sampling_region=$3

    bedtools intersect -a "$sampling_region" -b "$genes_complement" \
    | awk -F'\t' -v len="$len" '{ if (($3-$2) > len) print $0"\t"$3-$2}' \
    | sort -k1,1 -k2,2n -k3,3n \
    > "$complement_sampling_region" 
    
}

###########################
# Determines percentage of ambiguous nucleotides on sequences (used to filter out if >5% within get_sequence())
ambiguity_percent() {
    
    seq=$1

    amb_nucleotides=$( echo "$seq" | grep -o '[RYMKSWBDHVNrymkswbdhvn]' | wc -l )
    total_nucleotides=$( echo -n "$seq" | wc -c )
    ambiguity_percent=$(echo "scale=0; ($amb_nucleotides * 100) / $total_nucleotides" | bc)

    echo "$ambiguity_percent"

}

###########################
# Takes a random sample within a region:  
# selects a random start and end coordinate for the negative control from a specific region (closest to position=(exon+d))
sample_random() {
    chrRandom="$1"
    startRandom="$2"
    endRandom="$3"
    lenRandom="$4"

    # Generate random start coordinate within the specified range
    rand_start=$(( RANDOM % (endRandom - startRandom - lenRandom + 2) + startRandom ))      # 2 accounts for start being inclusive 

    # Calculate end coordinate based on start and region length
    rand_end=$(( rand_start + lenRandom - 1 ))

    # Print sampled region
    echo -e "$chrRandom\t$rand_start\t$rand_end" > "$random_output" 
}

###########################
# Extracts the sequence from the genome sequnece file, 
# calculates the distance of the negative control to the corresponding gene, 
# outputs the values on the final file 
get_sequence(){

    chrSeq=$1
    startSeq=$2
    endSeq=$3
    count=$4
    output_file=$5
    
    # Calculate distance to gene border 
    if [ "$flow" == 'downstream' ]; then distance_gene=$(( startClosest - end_gene )); 
    elif [ "$flow" == 'upstream' ]; then distance_gene=$(( start_gene - endClosest )); 
    fi

    # Extract sequence from genome and reject if it has an ambiguity percent > 5
    chromo=$( echo "$chrSeq" | tr -d "chr" )    
    seq=$( grep -w "chromosome $chromo" "$genome_seq" | cut -f 2 | cut -c"$startSeq"-"$endSeq" )
        
    if [ -n "$seq" ] && [[ $( ambiguity_percent "$seq" ) -lt 5 ]]; then                     

        countId=$(( count + $(wc -l < "$output_file") ))   

        echo "RNA$countId,No,$chrSeq,$startSeq,$endSeq,$seq,$distance_gene" >> "$output_file"

    fi

}

###########################
# Organize auxiliary functions and defines the positions=(exon+d)
#### For sequences with more than 1 exon #### 
multiple_exons_negative_control() {

    flow=$1
    sampling_regions=$2

    len_region=$(( lenA + lenB )) 

    for d in "${distances_to_seq[@]}"; do    
    
        # Downstream
        if [ "$flow" == 'downstream' ]; then
            
            initial_end=$(( end_gene + d ))
            initial_start=$(( initial_end - len_region ))
        
        # Upstream
        elif [ "$flow" == 'upstream' ]; then 

            initial_start=$(( start_gene - d ))
            if [ "$initial_start" -lt 0 ]; then initial_start=0; fi             # To stay within the chromosome region (downstream is limitated by the intersect of the complement)

            initial_end=$(( initial_start + len_region ))
        
        fi

    if [ "$initial_end" -gt 0 ] && [ "$initial_start" -gt 0 ]; then                         # Within chromosome limits  
            echo -e "$chr\t$initial_start\t$initial_end" > "$closest_input" 

            # Find closest region of genes complement 
            bedtools closest -a "$closest_input" -b "$sampling_regions" -D ref -t first  > "$closest_output"
        
            # Find random coordinates within the closest region
            IFS=$'\t' read -r _ _ _ chrClosest startClosest endClosest _ < "$closest_output" 
        
            #### EXON A ####
            sample_random "$chrClosest" "$startClosest" "$endClosest" "$lenA" 
            
            # Redefine variables with random coordinates 
            IFS=$'\t' read -r chrClosest startClosest endClosest < "$random_output" 

            get_sequence "$chrClosest" "$startClosest" "$endClosest" "$countA" "$negative_control_A"


            #### EXON B ####
            sample_random "$chrClosest" "$startClosest" "$endClosest" "$lenB" 
            
            # Redefine variables with random coordinates 
            IFS=$'\t' read -r chrClosest startClosest endClosest < "$random_output" 

            get_sequence "$chrClosest" "$startClosest" "$endClosest" "$countB" "$negative_control_B"
            
           
        fi
        
    done

    #rm -rf "$closest_input"
    #rm -rf "$closest_output"
}
    
#### For sequences with single exon (short-ncRNA) #### 
single_exon_negative_control() {

    flow=$1
    sampling_regions=$2
   
    for d in "${distances_to_seq[@]}"; do                                                       # Defined on line 15
    
        # Downstream
        if [ "$flow" == 'downstream' ]; then
            
            initial_end=$(( end_gene + d ))
            initial_start=$(( initial_end - len ))
        
        # Upstream
        elif [ "$flow" == 'upstream' ]; then 

            initial_start=$(( start_gene - d ))
            if [ "$initial_start" -lt 0 ]; then initial_start=0; fi 

            initial_end=$(( initial_start + len ))
        
        fi

        if [ "$initial_end" -gt 0 ] && [ "$initial_start" -gt 0 ]; then                         # Within chromosome limits  
            echo -e "$chr\t$initial_start\t$initial_end" > "$closest_input" 

            # Find closest region of genes complement 
            bedtools closest -a "$closest_input" -b "$sampling_regions" -D ref -t first  > "$closest_output"
        
            # Remove region which has already been used to avoid overlapping 
            IFS=$'\t' read -r _ _ _ chrClosest startClosest endClosest _ < "$closest_output" 
        
            sample_random "$chrClosest" "$startClosest" "$endClosest" "$len" 

            # Redefine variables with random coordinates 
            IFS=$'\t' read -r chrClosest startClosest endClosest < "$random_output" 

            get_sequence "$chrClosest" "$startClosest" "$endClosest" "$countSingle" "$negative_control_single"

        fi
        
    done

    #rm -rf "$closest_input"
    #rm -rf "$closest_output"
}
    
###########################

###########################################################################################################################
#### EXTRACT NEGATIVE CONTROLS ####
###########################################################################################################################
 
# Define if it's a single or multiple exon dataset counting the fields of initial data file (4: short-ncRNA; 5: lncRNA and protein-codingRNA)
num_fields=$( awk -F',' '{print NF}' "$initial_data" | sort -nu | tail -n 1 )                      

#### Single exon sequences ####
if [[ "$num_fields" -eq 4 ]]; then                                                                      

    if [ ! -s "$single_negative_coords" ]; then 

        echo "ID,Functional,Chromosome,Start,End,Sequence,Distance" > "$negative_control_single"
        countSingle=$(( $(wc -l < "$initial_data") - 1 ))                                           # Starts counting from last corresponding functional seq  
        
        tail -n +2 "$initial_data"  | while IFS=, read -r chr start_gene end_gene len; do   

          
            # Downstream 
            end_sampling_region=$(( end_gene + 5000000 ))
            echo -e "$chr\t$end_gene\t$end_sampling_region" > "$sampling_region"

            intersect "$sampling_region" "$len" "$complement_sampling_region"
           
            single_exon_negative_control 'downstream' "$complement_sampling_region"
            
            # Upstream 
            start_sampling_region=$(( start_gene - 5000000 ))
            if [ "$start_sampling_region" -lt 0 ]; then start_sampling_region=0; fi 

            echo -e "$chr\t$start_sampling_region\t$start_gene" > "$sampling_region"

            intersect "$sampling_region" "$len" "$complement_sampling_region"
           
            single_exon_negative_control 'upstream' "$complement_sampling_region" 

        done 

    fi

fi


# Multiple exons sequences #
if [[ "$num_fields" -eq 5 ]]; then 

    if [ ! -s "$first_negative_coords" ] || [ ! -s "$last_negative_coords" ]; then 

        echo "ID,Functional,Chromosome,Start,End,Sequence,Distance_to_functional" > "$negative_control_A"
        countA=$(( $(wc -l < "$initial_data") ))

        echo "ID,Functional,Chromosome,Start,End,Sequence,Distance_to_functional" > "$negative_control_B"
        countB=$(( $(wc -l < "$initial_data") ))

        tail -n +2 "$initial_data"  | while IFS=, read -r chr start_gene end_gene lenA lenB; do  

            len_region=$(( lenA+lenB ))
            # Temporary files
            sampling_region="$file_name"-sampling-region
            complement_sampling_region="$file_name"-complement-sampling-region

            # Downstream 
            end_sampling_region=$(( end_gene + 5000000 ))
            echo -e "$chr\t$end_gene\t$end_sampling_region" > "$sampling_region"

            intersect "$sampling_region" "$len_region" "$complement_sampling_region"
           
            multiple_exons_negative_control 'downstream' "$complement_sampling_region"
            
            # Upstream 
            start_sampling_region=$(( start_gene - 5000000 ))
            if [ "$start_sampling_region" -lt 0 ]; then start_sampling_region=0; fi 

            echo -e "$chr\t$start_sampling_region\t$start_gene" > "$sampling_region"

            intersect "$sampling_region" "$len_region" "$complement_sampling_region"
           
            multiple_exons_negative_control 'upstream' "$complement_sampling_region" 

        done 

    fi

fi



# Remove duplicates, sort sequences and rename (RNAid) 

format_output(){

    local file_with_duplicates=$1
    local output_file=$2

    local file_id="$(basename "${file_with_duplicates%.*}" | sed 's/-dataset//')" 
    local file_no_duplicates_unsorted="$file_id"-unsorted

    

    awk -F',' 'BEGIN { OFS="\t" }
            {
                key = $2 FS $3 FS $4 FS $5 FS $6
                count[key]++
                if (count[key] <= 1)
                    print $0
            }' "$file_with_duplicates" > "$file_no_duplicates_unsorted"


    awk -F ',' 'NR > 1 {print $2 "," $3 "," $4 "," $5 "," $6}' "$file_no_duplicates_unsorted" | sort -t ',' -k2,2 -k3,3n -k4,4n > "$file_id"-sorted-columns
    
    numb_seqs=$(wc -l < "$file_id"-sorted-columns)
    awk -F ',' 'NR > 1 {print $1}' "$file_with_duplicates" | head -n "$numb_seqs" > "$file_id"-id-column

    (echo "ID,Functional,Chromosome,Start,End,Sequence"; paste -d ',' "$file_id"-id-column "$file_id"-sorted-columns) > "$output_file"

}