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
#        $3 is the file which contains the sequence coordinates sorted and bed formatted
#        $4 is the location of the folder with the local databases/datasets or version specific executables
#       
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.

############################################################################################################################

initial_data=$1
initial_fasta=$2

first_rna_id=$(awk -F',' 'NR==2 {print $1}' "$initial_data" | tr -d 'RNA')
last_rna_id=$(awk -F',' 'END {print $1}' "$initial_data" | tr -d 'RNA') 

mkdir -p data/repeats
output_directory=data/repeats
output_file_copy="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')copy-number.csv"                  # Define name and directory for output file
output_file_distance="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dfam-distance.csv" 
 
echo "Genome_copy_number, T2T_copy_number" > "$output_file_copy"
echo "Chromosome,Start,End,Dfam_min,Dfam_sum" > "$output_file_distance"

##### Variables

mmseqs_exe=mmseqs
human_genome_mmseqs=data/raw/mmseqs/human_genome
T2T_genome_mmseqs=data/raw/mmseqs/T2T/human_genome         

bedtools_exe=bedtools
dfam_hits=data/raw/dfam-hg38-sorted.bed

############################################################################################################################

# Calculate genomic copy number (identify number of matches within human genome)

############################################################################################################################

## Run MMseqs2

# Hits agains hg38
mmseqs_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')mmseqs-output"

if [ ! -e "$mmseqs_output" ]; then
    $mmseqs_exe easy-search \
    "$initial_fasta" "$human_genome_mmseqs" "$mmseqs_output" "$output_directory/tmp" \
    --max-seqs 20000000000000000 --search-type 3 --format-output query,target,evalue \
    >/dev/null 2>>errors.log                              # max-seq not to be limited by default 600; search-type 3 for nucleotides; format-output to extract only required values 
fi

# Hits agains T2T
mmseqs_output_T2T="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')mmseqs-output-T2T"

if [ ! -e "$mmseqs_output_T2T" ]; then
    $mmseqs_exe easy-search \
    "$initial_fasta" "$T2T_genome_mmseqs" "$mmseqs_output_T2T" "$output_directory/tmp-t2t" \
    --max-seqs 20000000000000000 --search-type 3 --format-output query,target,evalue \
    >/dev/null 2>>errors.log
fi


## Process mmseqs results in order 

var=$first_rna_id
last_seq=$last_rna_id

while [ $var -lt $last_seq ]; do
    
    total_mmseqs=$( grep -w "RNA$var" $mmseqs_output | wc -l) 
    total_mmseqs_T2T=$( grep -w "RNA$var" $mmseqs_output_T2T | wc -l) 

    if [ -z "$total_mmseqs" ]; then total_mmseqs='NA'; fi 
    if [ -z "$total_mmseqs_T2T" ]; then total_mmseqs_T2T='NA'; fi

    echo "$total_mmseqs, $total_mmseqs_T2T" >> "$output_file_copy"

    (( var++ ))
   
done


############################################################################################################################

# Distance to nearest transposable element (Dfam)

############################################################################################################################

######## Format data for bedtools and temporary files

initial_data_temporary="${output_directory}/$(basename "${initial_data%.*}" | sed 's/-dataset//').bed" 
initial_data_sorted="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')sorted.bed"

awk -F',' 'NR > 1 {print $3, $4, $5}' "$initial_data" > "$initial_data_temporary"
sort -k1,1V -k2,2n "$initial_data_temporary" | tr ' ' '\t' > "$initial_data_sorted" 

downstream_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dfam-downstream.bed" 
upstream_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dfam-upstream.bed"
combined_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dfam-combined.bed"

######## Upstrean hits
$bedtools_exe closest -a "$initial_data_sorted" -b "$dfam_hits" -io -D ref -iu > "$downstream_output" 2>>errors.log

######## Downstream hits
$bedtools_exe closest -a "$initial_data_sorted" -b "$dfam_hits" -io -D ref -id > "$upstream_output" 2>>errors.log

paste <( cut -f 1,2,3,7 "$downstream_output" ) <( cut -f 7 "$upstream_output" ) --delimiters '\t' > "$combined_output"


######## Calculate sum of upstream and downstream, and combine into one file
while IFS=$'\t' read -r chr start end downstream upstream; do

    chr_numb=$(      echo "$chr" | cut -c 4- )
    upstream=$( echo "$upstream" | tr -d '-' )

    [ -z "$downstream" ] && downstream=0                                                # Implies that no region up or downstream, rather than data missing.
    [ -z "$upstream" ] && upstream=0 

    sum=$(( $downstream+$upstream ))

    if [[ "$upstream" -lt "$downstream" ]]; then                                        # Taking into account forward and reverse strand
                                               
        echo "$chr,$start,$end,$upstream,$sum" >> "$output_file_distance"

    else

        echo "$chr,$start,$end,$downstream,$sum" >> "$output_file_distance"   

    fi

done < "$combined_output"

#### Remove unnessesary files
rm -rf "$mmseqs_output"
rm -rf "$mmseqs_output_T2T"
rm -rf "$downstream_output"
rm -rf "$upstream_output"
rm -rf "$combined_output"
rm -rf "$initial_data_temporary"
rm -rf "$initial_data_sorted"

rm -rf "$output_directory"/tmp
rm -rf "$output_directory"/tmp-t2t



