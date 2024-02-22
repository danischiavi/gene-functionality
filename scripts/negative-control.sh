
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

file_name=data/$(basename "${initial_data%.*}" | sed 's/-coords-negative-control//')                 

# Temporary files
first_negative_coords="$file_name"-first-negative-coords.csv
last_negative_coords="$file_name"-last-negative-coords.csv
single_negative_coords="$file_name"-negative-coords.csv

# Final file
first_negative_control="$file_name"-first-negative-control-dataset.csv
last_negative_control="$file_name"-last-negative-control-dataset.csv
single_negative_control="$file_name"-negative-control-dataset.csv

genome_csv=data/raw/GRCh38.p14_genome.csv
gencode_bed=data/raw/gencode-ncrna-annotation.bed 
uniprot_bed=data/raw/unipAliSwissprot.bed

bedtools_exe=bedtools

###########################################################################################################################

#### Function declaration ####

#################################################################################################
# To extract coordinates and sequences of the negative control given a distance from the sequence
#################################################################################################
ambiguity_percent() {
    seq=$1

    amb_nucleotides=$( echo "$seq" | grep -o '[RYMKSWBDHVNrymkswbdhvn]' | wc -l )
    total_nucleotides=$( echo -n "$seq" | wc -c )
    ambiguity_percent=$(( (amb_nucleotides * 100) / total_nucleotides ))

    echo "$ambiguity_percent"

}

#### For sequences with more than 1 exon #### 
multiple_exons_negative_control() {

    local distance_to_seq=$1
    
    chromo=$( echo $chr | tr -d "chr" )                                                                                                    # Variables extract from initial_data file
    distance_between_exons=$(( "$end_seq" - "$start_seq" - "$first_len" - "$last_len" ))                                                    


    # UPSTREAM
    last_exon_left_end=$((     "$start_seq" - "$distance_to_seq" ))
    last_exon_left_start=$(( "$last_exon_left_end" - "$last_len" ))                                                                         # length based on exons's length
    
    first_exon_left_end=$(( "$last_exon_left_start" - "$distance_between_exons" ))                              
    first_exon_left_start=$(( "first_exon_left_end" - "$first_len" ))


    if [[ "$last_left_end" -gt 0 ]] || [[ "$last_left_start" -gt 0 ]]; then                                                                 # Filter out if negative coord
    
        last_left_sequence=$( grep -w "chromosome $chromo" "$genome_csv" | cut -f 2 | cut -c$last_exon_left_start-$last_exon_left_end )     # Extract seq from genome file                                                                                      # Blank to filter out if negative coord

        if [ ! -z "$last_left_sequence" ] && [[ $( ambiguity_percent "$last_left_sequence" ) -lt 5 ]]; then                                 # If no sequence extracted or if the ambiguous nucleotides are more than 5%, then remove

            echo "$chr,$last_exon_left_start,$last_exon_left_end,$last_left_sequence,$distance_to_seq" >> "$last_negative_coords"
        
        fi

    fi
   

    if [[ "$first_left_end" -gt 0 ]] || [[ "$first_left_start" -gt 0 ]]; then 
    
        first_left_sequence=$( grep -w "chromosome $chromo" "$genome_csv" | cut -f 2 | cut -c$first_exon_left_start-$first_exon_left_end )  
 
        if [ ! -z "$first_left_sequence" ] && [[ $( ambiguity_percent "$last_left_sequence" ) -lt 5 ]]; then  
            
            echo "$chr,$first_exon_left_start,$first_exon_left_end,$first_left_sequence,$distance_to_seq" >> "$first_negative_coords"
    
        fi

    fi
   

    ## DOWNSTREAM
    first_exon_right_start=$(( "$end_seq" + "$distance_to_seq" ))
    first_exon_right_end=$(( "$first_exon_right_start" + "$first_len" ))                                

    last_exon_right_start=$(( "$first_exon_right_end" + "$distance_between_exons" ))
    last_exon_right_end=$((  "$last_exon_right_start" + "$last_len" ))


    if [[ "$first_exon_right_start" -gt 0 ]] || [[ "$first_exon_right_end" -gt 0 ]]; then 
  
        first_right_sequence=$( grep -w "chromosome $chromo" $genome_csv | cut -f 2 | cut -c$first_exon_right_start-$first_exon_right_end ) 
    
        if [ ! -z "$first_right_sequence" ] && [[ $( ambiguity_percent "$last_left_sequence" ) -lt 5 ]]; then                                                                                       # If no sequence extracted, then remove.
            
            echo "$chr,$first_exon_right_start,$first_exon_right_end,$first_right_sequence,$distance_to_seq" >> "$first_negative_coords"
        
        fi
    
    fi
    

    if [[ "$last_exon_right_start" -gt 0 ]] || [[ "$last_exon_right_end" -gt 0 ]]; then 
  
        last_right_sequence=$( grep -w "chromosome $chromo" $genome_csv | cut -f 2 | cut -c$last_exon_right_start-$last_exon_right_end ) 
    
        if [ ! -z "$last_right_sequence" ] && [[ $( ambiguity_percent "$last_left_sequence" ) -lt 5 ]]; then                                                                                       # If no sequence extracted, then remove.
            
            echo "$chr,$last_exon_right_start,$last_exon_right_end,$last_right_sequence,$distance_to_seq" >> "$last_negative_coords"
        
        fi
    
    fi

}


#### For sequences with single exon (short-ncRNA) ####
single_exon_negative_control() {

    local distance_to_seq=$1
    
    chromo=$( echo $chr | tr -d "chr" )


    # UPSTREAM
    left_end=$(( "$start_seq" - "$distance_to_seq" ))
    left_start=$(( "$left_end" - "$len" ))                                                              
    
    if [[ "$left_end" -gt 0 ]] || [[ "$left_start" -gt 0 ]]; then 
    
        left_sequence=$( grep -w "chromosome $chromo" "$genome_csv" | cut -f 2 | cut -c$left_start-$left_end )
   
        if [ ! -z "$left_sequence"] && [[ $( ambiguity_percent "$last_left_sequence" ) -lt 5 ]]; then                                           # If no sequence extracted, then remove.

            echo "$chr,$left_start,$left_end,$left_sequence,$distance_to_seq" >> "$single_negative_coords"
        
        fi
    
    fi


    # DOWNSTREAM                       
    right_start=$(( "$end_seq" + "$distance_to_seq" ))
    right_end=$(( "right_start" + "$len" ))                                
                              

    if [[ "$right_start" -gt 0 ]] || [[ "$right_end" -gt 0 ]]; then 

        right_sequence=$( grep -w "chromosome $chromo" $genome_csv | cut -f 2 | cut -c$right_start-$right_end ) 

        if [ ! -z "$right_sequence" ] && [[ $( ambiguity_percent "$last_left_sequence" ) -lt 5 ]]; then   

            echo "$chr,$right_start,$right_end,$right_sequence,$distance_to_seq" >> "$single_negative_coords"

        fi
    
    fi

}



################################################################
# To filter negative control sequences using UniProt and GENCODE
################################################################

filter_out_functional(){
    
    local negative_coords=$1
    local exon=$2
    local negative_control=$3

    ## Temporary files
    bed_coords="$file_name"-"$exon"-coordinates.bed
    bed_intersect_uniprot="$file_name"-"$exon"-intersect-uniprot.bed
    bed_intersect_gencode="$file_name"-"$exon"-intersect-gencode.bed
    intersect_uniprot="$file_name"-"$exon"-intersect-uniprot
    intersect_gencode="$file_name"-"$exon"-intersect-gencode

    ######## Reformat data for bedtools
    grep -v "Start" "$negative_coords" | awk -F',' '{print $1 "\t" $2 "\t" $3}' > "$bed_coords"

    ######## Uniprot filtering
    "$bedtools_exe" intersect -a "$uniprot_bed" -b "$bed_coords" > "$bed_intersect_uniprot"

    ######## Extract chromosome coordinates of overlapping negative control sequences
    cut -f 1-3 "$bed_intersect_uniprot" | tr '\t' ',' > "$intersect_uniprot"
    #left_count=$( cat "$intersect" | wc -l )

    ######## GENCODE filtering
    "$bedtools_exe" intersect -a "$gencode_bed" -b "$bed_coords" > "$bed_intersect_gencode"

    ######## Extract chromosome coordinates of overlapping negative control sequences
    cut -f 1-3 "$bed_intersect_gencode" | tr '\t' ',' > "$intersect_gencode"

    ######## Remove negative control sequences that overlap with known functional sequences
    while IFS=',' read -r chr start end seq distance; do

        coordinates="${chr},${start},${end}"

        if grep -qF "$coordinates" "$intersect_uniprot"; then
            :
        elif grep -qF "$coordinates" "$intersect_gencode"; then
            :
        else
            (( count++ ))
            echo "RNA$count,No,$chr,$start,$end,$seq,$distance" >> "$negative_control"
        fi

    done < "$negative_coords"

    #rm -rf "$bed_coords"
    #rm -rf "$bed_intersect_uniprot"
    #rm -rf "$bed_intersect_gencode"
    #rm -rf "$intersect_uniprot"
    #rm -rf "$intersect_gencode"

}



###########################################################################################################################

# Generate sequences upstream and downstream of functional gene and filter out 

###########################################################################################################################

distances_to_seq=("1000" "10000" "100000" "1000000" "5000000")                             # Distances choosen for negative control (see manuscript for justification)

num_fields=$( awk -F',' '{print NF}' "$initial_data" | sort -nu | tail -n 1 )              # To define if it's single or multiple exon sequences


# Single exon sequences # 
if [[ "$num_fields" -eq 4 ]]; then                                                                      

    if [ ! -s "$single_negative_coords" ]; then 

        tail -n +2 "$initial_data"  | while IFS=, read chr start_seq end_seq len; do   

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

        tail -n +2 "$initial_data"  | while IFS=, read chr start_seq end_seq first_len last_len; do  

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
















######## Generate negative control FASTA file
#grep -v "Start" ./data/negative-control-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > ./data/negative-control-seq.fa

###########################################################################################################################

#rm -rf "$negative_coords"




