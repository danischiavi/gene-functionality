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
 
echo Genome_copy_number, mmseqs_copy_number, nhmmer_copy_number, T2T_copy_number > $output_file_copy
echo Chromosome,Start,End,Dfam_min,Dfam_sum > $output_file_distance


##### Variables
blastn_exe=blastn
human_genome_blastn=data/raw/blastn/human_genome

mmseqs_exe=mmseqs
human_genome_mmseqs=data/raw/mmseqs/human_genome

nhmmer_exe=nhmmer
human_genome_nhmmer=data/raw/GRCh38_p14_genomic.fna

T2T_genome_mmseqs=data/raw/mmseqs/T2T/human_genome         

bedtools_exe=bedtools
dfam_hits=data/raw/dfam-hg38-nrph.bed

############################################################################################################################

# Calculate genomic copy number (identify number of matches within human genome)

############################################################################################################################

## Run blastn and MMseqs2

blastn_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')blastn-output.csv" 
$blastn_exe -query $initial_fasta -db $human_genome_blastn -evalue 0.01 -out $blastn_output -outfmt "10 qaccver saccver pident" >/dev/null 2>>errors.log

mmseqs_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')mmseqs-output"
$mmseqs_exe easy-search "$initial_fasta" "$human_genome_mmseqs" "$mmseqs_output" ""$output_directory"/tmp" --search-type 3 --format-output query,target,evalue >/dev/null 2>>errors.log

mmseqs_output_T2T="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')mmseqs-output-T2T"
$mmseqs_exe easy-search "$initial_fasta" "$T2T_genome_mmseqs" "$mmseqs_output_T2T" ""$output_directory"/tmp-t2t" --search-type 3 --format-output query,target,evalue >/dev/null 2>>errors.log

## Run nhmmer and process Blastn and mmseqs results in order 

nhmmer_input="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')nhmmer-input.fasta"
nhmmer_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')nhmmer-output"

var=$first_rna_id
last_seq=$last_rna_id

while [ $var -lt $last_seq ]; do

    echo "$( grep -w -A 1 ">RNA$var" $initial_fasta )" > $nhmmer_input

    $nhmmer_exe --tblout $nhmmer_output -E 0.01 --noali $nhmmer_input $human_genome_nhmmer > /dev/null 2>&1 # -E <x> : report sequences <= this E-value threshold in output; don't output alignments
    nhmmer_hits=$( grep -v "#" $nhmmer_output | wc -l )

    
    #minimap2 data/raw/GRCh38_p14_genomic.fna data/nhmmer-input.fasta >> data/minimap-output
    
    ######## Process each RNA in the file (by ID) as the counter increases
    total_blastn=$( grep -w "RNA$var" $blastn_output | wc -l)  
    total_mmseqs=$( grep -w "RNA$var" $mmseqs_output | wc -l) 
    total_mmseqs_T2T=$( grep -w "RNA$var" $mmseqs_output_T2T | wc -l) 

    if [ -z "$total_blastn" ]; then total_blastn='NA'; fi 
    if [ -z "$total_mmseqs" ]; then total_mmseqs='NA'; fi 
    if [ -z "$nhmmer_hits" ]; then nhmmer_hits='NA'; fi
    if [ -z "$total_mmseqs_T2T" ]; then total_mmseqs_T2T='NA'; fi

    echo "$total_blastn, $total_mmseqs, $nhmmer_hits", "$total_mmseqs_T2T" >> $output_file_copy

    (( var++ ))
   
done

## Telomere-to-telomere assembly with nhmmer

#var=$first_rna_id
#last_seq=$last_rna_id

#while [ $var -le $last_seq ]
#do
#    echo "$( grep -w -A 1 ">RNA$var" $initial_fasta )" > data/nhmmer-input.fasta
#    $nhmmer_exe --tblout data/nhmmer-output -E 0.01 --noali data/nhmmer-input.fasta $T2T_genome > /dev/null 2>&1 # -E <x> : report sequences <= this E-value threshold in output; don't output alignments
#    nhmmer_hits_T2T=$( grep -v "#" data/nhmmer-output | wc -l )
    
#    if [ -z "$nhmmer_hits_T2T" ]; then nhmmer_hits_T2T='NA'; fi

#    echo "$nhmmer_hits_T2T" >> data/T2T-copy-number-$name.csv

#    (( var++ ))

#done



############################################################################################################################

# Distance to nearest transposable element (Dfam)

############################################################################################################################

######## Format data for bedtools and temporary files

initial_data_temporary="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//').bed" 
initial_data_sorted="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')sorted.bed"

awk -F',' 'NR > 1 {print $3, $4, $5}' "$initial_data" | tr ',' '\t' > $initial_data_temporary
sort -k1,1 -k2,2n -o $initial_data_sorted $initial_data_temporary


downstream_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dfam-downstream.bed" 
upstream_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dfam-upstream.bed"
combined_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dfam-combined.bed"

######## Upstrean hits
$bedtools_exe closest -a $initial_data_sorted -b $dfam_hits -io -D ref -iu > $downstream_output 2>>errors.log

######## Downstream hits
$bedtools_exe closest -a $initial_data_sorted -b $dfam_hits -io -D ref -id > $upstream_output 2>>errors.log

paste <( cut -f 1,2,3,7 $downstream_output ) <( cut -f 7 $upstream_output ) --delimiters '\t' > $combined_output


######## Calculate sum of upstream and downstream, and combine into one file
while IFS=$'\t' read -r chr start end downstream upstream; do

    chr_numb=$(      echo "$chr" | cut -c 4- )
    upstream=$( echo "$upstream" | tr -d '-' )

    [ -z "$downstream" ] && downstream=0                                                # Implies that no region up or downstream, rather than data missing.
    [ -z "$upstream" ] && upstream=0 

    sum=$(( $downstream+$upstream ))

    if [[ "$upstream" -lt "$downstream" ]]; then                                        # Taking into account forward and reverse strand
                                               
        echo "$chr,$start,$end,$upstream,$sum" >> $output_file_distance

    else

        echo "$chr,$start,$end,$downstream,$sum" >> $output_file_distance   

    fi

done < $combined_output

#### Remove unnessesary files
rm -rf $downstream_output
rm -rf $upstream_output
rm -rf $combined_output
rm -rf $initial_data_temporary
#rm -rf $initial_data_sorted

rm -rf data/nhmmer-input.fasta
rm -rf data/nhmmer-output
# rm -rf ""$output_directory"/tmp"




