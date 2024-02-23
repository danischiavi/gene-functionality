#!/bin/bash
#
# Script Name: genomic-repeat-associated-features.sh
#
#
# Description: This script calculates: 
#                       - minimum and average distance to repeat 
#                       - closest non overlapping repetitive element 
#                       - genomic copy number  
# for protein-coding exons, lncRNA exons, short ncRNAs and negative control sequences. 
#
# Input: $1 is the dataset 
#        $2 is the fasta file for the sequences
#        $3 
#        
#       
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.

############################################################################################################################

#### Files and directories #### 
initial_data=$1
initial_fasta=$2
human_genome_mmseqs=$3
dfam_hits=$4

first_rna_id=$(awk -F',' 'NR==2 {print $1}' "$initial_data" | tr -d 'RNA')
last_rna_id=$(awk -F',' 'END {print $1}' "$initial_data" | tr -d 'RNA') 

output_directory=data/repeats
mkdir -p "$output_directory"

file_name="$output_directory"/$(basename "${initial_data%.*}" | sed 's/-dataset//')

# Temporary files
output_file_copy="$file_name"-copy-number.csv                                        
output_file_distance_not_matching="$file_name"-dfam-distance-not-matching.csv
output_file_distance_temp="$file_name"-dfam-distance-temp.csv
output_file_distance="$file_name"-dfam-distance.csv

# Final output file 
output_file_repeats="$file_name"-repeats.csv 

##### Variables #####
mmseqs_exe=mmseqs        
bedtools_exe=bedtools

############################################################################################################################

# Calculate genomic copy number (identify number of matches within human genome)

############################################################################################################################

if [ ! -s "$output_file_copy" ]; then

    #### Run MMseqs2 ####

    # Hits agains hg38
    mmseqs_output="$file_name"-mmseqs-output

    if [ ! -e "$mmseqs_output" ]; then

        "$mmseqs_exe" easy-search "$initial_fasta" "$human_genome_mmseqs" "$mmseqs_output" "$file_name-tmp" --max-seqs 20000000000000000 --search-type 3 --format-output query,target,evalue >/dev/null 2>>errors.log                              # max-seq not to be limited by default 600; search-type 3 for nucleotides; format-output to extract only required values 
    
    fi

    ## Process mmseqs results in order 

    echo "copy_number" > "$output_file_copy"

    var=$first_rna_id
    last_seq=$last_rna_id

    while [ $var -le $last_seq ]; do
    
        total_mmseqs=$( grep -w "RNA$var" $mmseqs_output | wc -l) 

        if [ -z "$total_mmseqs" ] || [[ "$total_mmseqs" == 0 ]]; then 
        
            total_mmseqs='NA'
        
        fi 
        
        echo "$total_mmseqs" >> "$output_file_copy"

        (( var++ ))
   
    done

fi

############################################################################################################################

# Distance to nearest transposable element (Dfam)

############################################################################################################################

if [ ! -s "$output_file_distance" ]; then

    echo "Dfam_min,Dfam_sum" > "$output_file_distance"

    ##### Format data for bedtools and temporary files ####
    # Bedtools requires sorted input which matches the database ($dfam_hits)

    initial_data_temporary="$file_name".bed 
    initial_data_sorted="$file_name"-sorted.bed

    awk -F',' 'NR > 1 {print $3, $4, $5}' "$initial_data" > "$initial_data_temporary"
    sort -k1,1V -k2,2n "$initial_data_temporary" | tr ' ' '\t' > "$initial_data_sorted" 

    downstream_output="$file_name"-dfam-downstream.bed
    upstream_output="$file_name"-dfam-upstream.bed
    combined_output="$file_name"-dfam-combined.bed

    # Upstrean hits
    $bedtools_exe closest -a "$initial_data_sorted" -b "$dfam_hits" -io -iu -D ref  > "$downstream_output" 2>>errors.log 
    # -io : ignore features in B that overlap A
    # -iu : ignore upstream Ignore features in B that are upstream of features in A. Required -D option  
    # -D: report the closest feature in B, and its distance to A as an extra column

    # Downstream hits
    $bedtools_exe closest -a "$initial_data_sorted" -b "$dfam_hits" -io -id -D ref  > "$upstream_output" 2>>errors.log

    # Combine outputs 
    paste <( cut -f 1,2,3,7 "$downstream_output" ) <( cut -f 7 "$upstream_output" ) --delimiters '\t' > "$combined_output"


    ## Calculate sum of upstream and downstream, and combine into one file
    while IFS=$'\t' read -r chr start end downstream upstream; do

        chr_numb=$(      echo "$chr" | cut -c 4- )
        upstream=$( echo "$upstream" | tr -d '-' )

        [ -z "$downstream" ] && downstream=0                                                # Implies that no region up or downstream, rather than data missing.
        [ -z "$upstream" ] && upstream=0 

        sum=$(( $downstream+$upstream ))

        if [[ "$upstream" -lt "$downstream" ]]; then                                        # Taking into account forward and reverse strand
                                               
            echo "$chr,$start,$end,$upstream,$sum" >> "$output_file_distance_not_matching"

        else

            echo "$chr,$start,$end,$downstream,$sum" >> "$output_file_distance_not_matching"   

        fi

    done < "$combined_output"

fi


#### Reformat intermediate output files ####
# Bedtools output has to be sorted again to match the initial data order

# Rscript: Matches start coordinate and chr, and returns a file with matching bedtools results and initial_data following initial_data order 

Rscript scripts/reformat_dfam.R "$initial_data" "$output_file_distance_not_matching" "$output_file_distance_temp" >/dev/null 2>>errors.log

cat "$output_file_distance_temp" | cut -d ',' -f 7,8 > "$output_file_distance"


## Join output files for better organization 
if [ ! -s "$output_file_repeats" ]; then

    paste -d',' "$output_file_copy" "$output_file_distance" > "$output_file_repeats"

fi


## Remove unnessesary files
rm -rf "$mmseqs_output"
rm -rf "$downstream_output"
rm -rf "$upstream_output"
rm -rf "$combined_output"
rm -rf "$initial_data_temporary"
rm -rf "$initial_data_sorted"
rm -rf "$output_file_distance_not_matching"
rm -rf "$output_file_distance_temp"


#rm -rf "$output_file_copy"
#rm -rf "$output_file_distance"


