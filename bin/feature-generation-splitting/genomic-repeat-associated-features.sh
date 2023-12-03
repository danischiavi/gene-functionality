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
sorted-coordinates.bed=$3
additional_folder=$4

####Checks 

## Local database for blastn

blast_total=$( ls $additional_folder/human_genome* | wc -l )
if [[ $blast_total -ge 6 ]]
then
    human_genome=$additional_folder/human_genome
else
    until [ -f $human_genome.nhr ] ; do read -p "Please enter custom name of human genome blastn database: " h ; human_genome=$additional_folder/$h ; echo ; done
fi

## blastn

if find $additional_folder/ -executable -type f | grep -q blastn
then
    blastn_exe=$additional_folder/blastn
elif command -v blastn &> /dev/null
then
    blastn_exe=blastn
else
    echo Please check that blastn is installed.
    exit 1
fi


#### File generation
echo Genome_copy_number,Genome_complete_match > copy-number.csv
echo Chromosome,Start,End,Dfam_min,Dfam_sum > dfam-distance.csv

############################################################################################################################

# Calculate genomic copy number

############################################################################################################################

######## Run blastn to identify number of matches within human genome
$blastn_exe -query $initial_fasta -db $human_genome -evalue 0.01 -out output.csv -outfmt "10 qaccver saccver pident" >/dev/null 2>>errors.log

######## Counter to process each RNA in order
count_file=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
total_rna=$( tail -1 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )

while [ $count_file -le $total_rna ]
do
    ######## Process each RNA in the file (by ID) as the counter increases
    grep -w "RNA$count_file" output.csv > results  
    total=$(cat results | wc -l)
    [ -z "$total" ] && echo NA >> copy-number.csv || echo $total >> copy-number.csv
    count_file=$(( $count_file + 1 ))
done

rm -rf results
rm -rf rna.fa
rm -rf output.csv
rm -rf sequences

############################################################################################################################

# Distance to nearest transposable element (Dfam)

############################################################################################################################

######## Upstrean hits
$bedtools_exe closest -a sorted-coordinates.bed -b $dfam_hits -io -D ref -iu > dfam_downstream.bed 2>>errors.log

######## Downstream hits
$bedtools_exe closest -a sorted-coordinates.bed -b $dfam_hits -io -D ref -id > dfam_upstream.bed 2>>errors.log

paste <( cut -f 1,2,3,7 dfam_downstream.bed ) <( cut -f 7 dfam_upstream.bed ) --delimiters '\t' > dfam_combined.bed


######## Calculate sum of upstream and downstream, and combine into one file
while IFS=$'\t' read -r chr start end downstream upstream
do
    chr_numb=$( echo "$chr" | cut -c 4- )
    upstream=$( echo "$upstream" | tr -d '-' )
    [ -z "$downstream" ] && downstream=0    # Implies that no region up or downstream, rather than data missing.
    [ -z "$upstream" ] && upstream=0 
    sum=$(( $downstream+$upstream ))
    if [[ "$upstream" -lt "$downstream" ]]    # Taking into account forward and reverse strand
    then
        echo $chr,$start,$end,$upstream,$sum >> dfam-distance.csv
    elif [[ "$downstream" -lt "$upstream" ]]
    then
        echo $chr,$start,$end,$downstream,$sum >> dfam-distance.csv
    else
        echo $chr,$start,$end,$downstream,$sum >> dfam-distance.csv #??? is it necessary? 
    fi

done < dfam_combined.bed

mv dfam_downstream.bed additional-output/
mv dfam_upstream.bed additional-output/
rm -rf dfam_combined.bed