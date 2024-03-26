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
output_file_distance_not_sorted="$file_name"-dfam-distance-not-sorted.csv
output_file_distance_tmp="$file_name"-dfam-distance-tmp.csv
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

    echo "repeat_distance" > "$output_file_distance"

    ## Temporary files ## 
    downstream_output="$file_name"-dfam-downstream.bed
    upstream_output="$file_name"-dfam-upstream.bed
    combined_output="$file_name"-dfam-combined.bed
    dfam_hits_tmp="$file_name"-dfam-hits-tmp
    initial_data_tmp="$file_name".bed 
    initial_data_sorted="$file_name"-sorted.bed

    ##### Format data for bedtools ####
    # Bedtools requires sorted input which matches the database ($dfam_hits)
    awk -F',' 'NR > 1 {print $3, $4, $5}' "$initial_data" > "$initial_data_tmp"
    sort -k1,1V -k2,2n "$initial_data_tmp" | tr ' ' '\t' > "$initial_data_sorted" 

    # Bedtools closest doesn't work if any chr is missing. This is unlikely on a big dataset, but to avoid errors: lines which corresponds to chr not present on the dataset are filtered out of the Dfam_hits database 
    awk -F'\t' 'NR==FNR{chrom[$1]; next} $1 in chrom' "$initial_data_sorted" "$dfam_hits" > "$dfam_hits_tmp"

    
    # Upstrean hits
    $bedtools_exe closest -a "$initial_data_sorted" -b "$dfam_hits_tmp" -io -iu -D ref  > "$downstream_output" 2>>errors.log 
    # -io : ignore features in B that overlap A
    # -iu : ignore upstream Ignore features in B that are upstream of features in A. Required -D option  
    # -D: report the closest feature in B, and its distance to A as an extra column

    # Downstream hits
    $bedtools_exe closest -a "$initial_data_sorted" -b "$dfam_hits_tmp" -io -id -D ref  > "$upstream_output" 2>>errors.log

    # Combine outputs 
    paste <( cut -f 1,2,3,7 "$downstream_output" ) <( cut -f 7 "$upstream_output" ) --delimiters '\t' > "$combined_output"


    ## Calculate sum of upstream and downstream, and combine into one file
    while IFS=$'\t' read -r chr start end downstream upstream; do

        chr_numb=$(      echo "$chr" | cut -c 4- )
        upstream=$( echo "$upstream" | tr -d '-' )

        [ -z "$downstream" ] && downstream=0                                                # Implies that no region up or downstream, rather than data missing.
        [ -z "$upstream" ] && upstream=0 

        sum=$(( $downstream+$upstream ))

        echo "$chr,$start,$end,$sum" >> "$output_file_distance_not_sorted"


    done < "$combined_output"

fi


#### Reformat intermediate output files ####
# Bedtools output has to be sorted again to match the initial data order

# Rscript: Matches start coordinate and chr, and returns a file with matching bedtools results and initial_data following initial_data order 

Rscript scripts/reformat_dfam.R "$initial_data" "$output_file_distance_not_sorted" "$output_file_distance_tmp" >/dev/null 2>>errors.log

awk -F',' '{print $(NF-1)","$NF}' "$output_file_distance_tmp" > "$output_file_distance"             # last 2 columns have the values (not same column number for negative and functional dataset due to distance label)

## Join output files for better organization 
if [ ! -s "$output_file_repeats" ]; then

    paste -d',' "$output_file_copy" "$output_file_distance" > "$output_file_repeats"

fi


## Remove unnessesary files
rm -rf "$mmseqs_output"
rm -rf "$downstream_output"
rm -rf "$upstream_output"
rm -rf "$combined_output"
rm -rf "$initial_data_tmp"
rm -rf "$initial_data_sorted"
rm -rf "$output_file_distance_not_sorted"
rm -rf "$output_file_distance_tmp"
rm -rf "$file_name-tmp"
rm -rf "$dfam_hits_tmp"
#rm -rf "$output_file_copy"
#rm -rf "$output_file_distance"


