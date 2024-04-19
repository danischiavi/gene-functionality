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

distances_to_seq=("100" "1000" "10000" "100000" "1000000" "5000000")                             # Distances choosen for negative control (see manuscript for justification)
distance_between_exons=100

file_name=data/datasets/$(basename "${initial_data%.*}" | sed 's/-coords-negative-control//')                 

# Temporary files
first_negative_coords="$file_name"-first-negative-coords.csv
last_negative_coords="$file_name"-last-negative-coords.csv
single_negative_coords="$file_name"-negative-coords.csv

# Output files
negative_control_single="$file_name"-negative-control-dataset.csv   # File for short-ncrna is the final negative control dataset
negative_control_first="$file_name"-first-negative-control.csv      # File for lncrna and protein will be rename after this script to associated with the corresponding exon  
negative_control_last="$file_name"-last-negative-control.csv

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

    # Temporary files                                                                                                                                     
    closest_input_unsorted="$file_name"-closest-input-unsorted
    closest_input="$file_name"-closest-input
    closest_output="$file_name"-closest-output
    closest_output_filtered="$file_name"-closest-output-filtered

    len_region=$(( first_len + distance_between_exons + last_len ))

    for d in "${distances_to_seq[@]}"; do                                                       # Defined on line 15
    
    # Upstream
        initial_end_up=$(( start_gene - d ))
        initial_start_up=$(( initial_end_up - len_region ))
        
        if [ "$initial_end_up" -gt 0 ] && [ "$initial_start_up" -gt 0 ]; then
            echo -e "$chr\t$initial_start_up\t$initial_end_up" >> "$closest_input_unsorted" 
        fi

    # Downstream
        initial_start_down=$(( end_gene + d ))
        initial_end_down=$(( initial_start_down + len_region ))

        if [ "$initial_end_down" -gt 0 ] && [ "$initial_start_down" -gt 0 ]; then
            echo -e "$chr\t$initial_start_down\t$initial_end_down" >> "$closest_input_unsorted" 
        fi

    done

    if [ -s "$closest_input_unsorted" ]; then 

        sort -k1,1 -k2,2n -k3,3n "$closest_input_unsorted" > "$closest_input" 

        # Find the closest region of the genes complement 
        bedtools closest -a "$closest_input" -b "$genes_complement" -D ref -t first  > "$closest_output"

        # Remove closest regions which are duplicated and don't have enough nucleotides to extract sequences for both exons 
        awk -v len_region="$len_region" -F'\t' 'BEGIN {OFS="\t"} {
            key = $4 "\t" $5 "\t" $6  
            count[key]++     
            if (count[key] <= 1 && $6 - $5 >= len_region) {  # Check if count is less than or equal to 1 and the region length condition is met
                print $0, $6 - $5  
            }
        }' "$closest_output" > "$closest_output_filtered"

        ## Format final output file
        while IFS=$'\t' read -r chr initial_start _ _ start end _ _; do

            chromo=$( echo "$chr" | tr -d "chr" )

            if [ "$start" -ge "$initial_start" ]; then               # if closest is downstream initial coords
        
                # First exon
                start_first="$start"
                end_first=$(( start + first_len ))

                # Last exon
                start_last=$(( end_first + distance_between_exons ))
                end_last=$(( start_last + last_len ))
        
            else                                                    # if closest is upstream initial coords
                                                                    
                # Last exon
                end_last="$end"
                start_last=$(( end_last - last_len ))

                # First exon
                end_first=$(( start_last - distance_between_exons )) 
                start_first=$(( end_first - first_len ))
        
            fi

            # Calculate distance to gene border 
            if [ "$start_first" -ge "$start_gene" ]; then           # if negative control is downstream gene

                distance_gene_first=$(( start_first - end_gene ))
                distance_gene_last=$(( start_last - end_gene )) 

            else                                                    # if negative control is upstream gene

                distance_gene_first=$(( start_gene - start_first ))
                distance_gene_last=$(( start_gene - start_last )) 
            
            fi

            # Extract sequence from genome and reject if it has an ambiguity percent > 5
            seq_first=$( grep -w "chromosome $chromo" "$genome_seq" | cut -f 2 | cut -c"$start_first"-"$end_first" )
        
            if [ -n "$seq_first" ] && [[ $( ambiguity_percent "$seq_first" ) -lt 5 ]]; then                     

                (( count_first++ ))
                echo "RNA$count_first,No,$chr,$start_first,$end_first,$seq_first,$distance_gene_first" >> "$negative_control_first"
        
            fi

            seq_last=$( grep -w "chromosome $chromo" "$genome_seq" | cut -f 2 | cut -c"$start_last"-"$end_last" )
        
            if [ -n "$seq_last" ] && [[ $( ambiguity_percent "$seq_last" ) -lt 5 ]]; then  

                (( count_last++ ))
                echo "RNA$count_last,No,$chr,$start_last,$end_last,$seq_last,$distance_gene_last" >> "$negative_control_last"
        
            fi

        done < "$closest_output_filtered"
    fi

    rm -rf "$closest_input_unsorted"
    rm -rf "$closest_input"
    rm -rf "$closest_output"
    rm -rf "$closest_output_filtered"
}
  

#### For sequences with single exon (short-ncRNA) 
single_exon_negative_control() {

    # Temporary files                                                                                                                                     
    closest_input_unsorted="$file_name"-closest-input-unsorted
    closest_input="$file_name"-closest-input
    closest_output="$file_name"-closest-output
    closest_output_filtered="$file_name"-closest-output-filtered

    len_region="$len"

    for d in "${distances_to_seq[@]}"; do                                                       # Defined on line 15
    
    # Upstream
        initial_end_up=$(( start_gene - d ))
        initial_start_up=$(( initial_end_up - len_region ))
        
        if [ "$initial_end_up" -gt 0 ] && [ "$initial_start_up" -gt 0 ]; then
            echo -e "$chr\t$initial_start_up\t$initial_end_up" >> "$closest_input_unsorted" 
        fi

    # Downstream
        initial_start_down=$(( end_gene + d ))
        initial_end_down=$(( initial_start_down + len_region ))

        if [ "$initial_end_down" -gt 0 ] && [ "$initial_start_down" -gt 0 ]; then
            echo -e "$chr\t$initial_start_down\t$initial_end_down" >> "$closest_input_unsorted" 
        fi

    done

    if [ -s "$closest_input_unsorted" ]; then 
    
        sort -k1,1 -k2,2n -k3,3n "$closest_input_unsorted" > "$closest_input" 

        # Find closest region of genes complement 
        bedtools closest -a "$closest_input" -b "$genes_complement" -D ref -t first  > "$closest_output"

        # Remove closest regions which are duplicated and don't have enough nucleotides to extract sequence
        awk -v len_region="$len_region" -F'\t' 'BEGIN {OFS="\t"} {
            key = $4 "\t" $5 "\t" $6  
            count[key]++    

            if (count[key] <= 1 && $6 - $5 >= len_region) {  # Check if count is less than or equal to 1 and the region length condition is met
                print $0, $6 - $5  
            }
        }' "$closest_output" > "$closest_output_filtered"

        ## Format final output file
        while IFS=$'\t' read -r chr initial_start _ _ start end _ _; do

            chromo=$( echo "$chr" | tr -d "chr" )

            if [ "$start" -ge "$initial_start" ]; then               # if closest is downstream initial coords
        
                start_single="$start"
                end_single=$(( start_single + len ))

        
            else                                                    # if closest is upstream initial coords
                                                        
                end_single="$end"
                start_single=$(( end_single - len ))
        
            fi

        
            # Calculate distance to gene border 

            if [ "$start_single" -ge "$start_gene" ]; then           # if negative control is downstream gene

                distance_gene=$(( start_single - end_gene )) 

            else                                                    # if negative control is upstream gene

                distance_gene=$(( start_gene - start_single ))
    
            fi

            # Extract sequence from genome and reject if it has an ambiguity percent > 5
        
            seq=$( grep -w "chromosome $chromo" "$genome_seq" | cut -f 2 | cut -c"$start_single"-"$end_single" )
        
            if [ -n "$seq" ] && [[ $( ambiguity_percent "$seq" ) -lt 5 ]]; then                     

                (( count_single++ ))
                echo "RNA$count_single,No,$chr,$start_single,$end_single,$seq,$distance_gene" >> "$negative_control_single"
        
            fi

        done < "$closest_output_filtered"
    fi

    rm -rf "$closest_input_unsorted"
    rm -rf "$closest_input"
    rm -rf "$closest_output"
    rm -rf "$closest_output_filtered"
}
    
###########################

###########################################################################################################################
#### EXTRACT NEGATIVE CONTROLS ####
###########################################################################################################################

num_fields=$( awk -F',' '{print NF}' "$initial_data" | sort -nu | tail -n 1 )                  # To define if it's single or multiple exon sequences

# Single exon sequences # 
if [[ "$num_fields" -eq 4 ]]; then                                                                      

    if [ ! -s "$single_negative_coords" ]; then 

        echo "ID,Functional,Chromosome,Start,End,Sequence,Distance" > "$negative_control_single"
        count_single=$(( $(wc -l < "$initial_data") - 1 ))                                     # Start counting from last functional seq  
        
        tail -n +2 "$initial_data"  | while IFS=, read -r chr start_gene end_gene len; do   

            single_exon_negative_control

        done 

    fi

fi


# Multiple exons sequences #
if [[ "$num_fields" -eq 5 ]]; then 

    if [ ! -s "$first_negative_coords" ] || [ ! -s "$last_negative_coords" ]; then 

        echo "ID,Functional,Chromosome,Start,End,Sequence,Distance_to_functional" > "$negative_control_first"
        count_first=$(( $(wc -l < "$initial_data") ))

        echo "ID,Functional,Chromosome,Start,End,Sequence,Distance_to_functional" > "$negative_control_last"
        count_last=$(( $(wc -l < "$initial_data") ))

        tail -n +2 "$initial_data"  | while IFS=, read -r chr start_gene end_gene first_len last_len; do  

            multiple_exons_negative_control

        done 

    fi

fi

