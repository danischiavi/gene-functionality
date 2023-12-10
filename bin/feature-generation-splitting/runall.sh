#!/bin/bash
#
# Script Name: runall-features.sh
#
#
# Description: This script runs all the scripts to calculates the chosen features for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences. 
#
# Input: $1 is the file identifier for the dataset and fasta file. Eg: 200702-functional-ncrna
#        $2 is the location of the folder with the local databases/datasets or version specific executables
# 
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.


############################################################################################################################

# Initialization 

############################################################################################################################

#set -xv
IFS=$'\n'
date=$(date +%y%m%d)

######## Function for calculation with float numbers

calc() { awk "BEGIN{print $*}"; }


######## Delete previous log files

rm -rf tabix.log
rm -rf errors.log


#### Global variables 

######## CSV and FASTA for initial data
initial_data=$1-dataset.csv 
initial_fasta=$1-seq.fa 

if [ ! -f $initial_data ] || [ ! -f $initial_fasta ]
then
    echo "Initial data files do not exist."
    exit 1
else
    :
fi


######## Folder with local databases and version specific executables
#additional_folder=$2  # is the best way to organize the project? 
#if [ ! -d $additional_folder/ ]
#then
#    echo "Folder with required local databases does not exist."
    #exit 1
#else
#    :
#fi

last_rna_id=$( tail -1 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
first_rna_id=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )

#### 
lncrna_fasta="$1"
short_ncrna_fasta="$2"
chr_coord_ncrna="$3"
ids_HGNC="$4"
genome_gff="$5"
genome_csv="$6"

# Input:
#        $1 is the downloaded fasta file from RNAcentral containing the lncRNAs. 
#        $2 is the downloaded fasta file from RNAcentral containing the short-ncRNAs.
#        $3 is the downloaded chromosome coordinates file from RNAcentral for all ncRNAs in the database. Eg: Homo_sapiens.GRCh38.bed
#        $4 is the text file containing RefSeq IDs from HGNC
#        $5 is the GFF file for the downloaded human genome Eg: GCF_000001405.39_GRCh38.p13_genomic.gff
#        $6 is the converted CSV file of the human genome. Eg: GRCh38.p13_genome.csv
##### 


######## Create folder for files generated

[ -d results/"$date"/ ] && rm -rf results/"$date"/* || mkdir -p results/"$date"/
echo



#### Reformat initial dataset
awk -F',' 'NR > 1 {print $3, $4, $5}' $initial_data > ./data/coordinates  # for mafFetch
awk -F',' 'NR > 1 {print $3":"$4"-"$5}' $initial_data > locations
awk -F',' 'NR > 1 {print $3"\t"$4"\t"$5}' "$initial_data" | sort -k1,1 -k2,2n > ./data/sorted-coordinates.bed # for bedtools
awk -F',' 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1}' "$initial_data" > ./data/hg38-coordinates.bed # for VFC 1kGP files


############################################################################################################################

# Calculate features

############################################################################################################################

./intrinsic-seq-features.sh 

./seq-conservation-features.sh 

./transcriptome-expression-features.sh 

./genomic-repeat-associated-features.sh 

./protein-and-rna-specific-features.sh

./population-variation-features.sh 

############################################################################################################################

# Combine all generated data into one file

############################################################################################################################

######## Update data parameter, in case script runs over multiple days
date=$(date +%y%m%d)
name=$( echo $initial_data | cut -d '-' -f 2,3 )

paste -d ',' $initial_data 1000g-freqsummary.csv rscape-dataset.csv interaction-intermediate.csv copy-number.csv tissue-rnaseq.csv primary-cell-rnaseq.csv GC.csv MFE-final.csv access.csv rnacode.csv ncrna-dfam-distance.csv conservation.csv CPC2.csv gerp.csv gnomad.csv > $date-$name-final-dataset.csv

######## Move previous files to folder, in case there was an issue during feature calculation
mv 1000g-freqsummary.csv additional-output/
mv rscape-dataset.csv additional-output/
mv interaction-intermediate.csv additional-output/
mv copy-number.csv additional-output/
mv tissue-rnaseq.csv additional-output/
mv primary-cell-rnaseq.csv additional-output/
mv GC.csv additional-output/
mv MFE-final.csv additional-output/
mv access.csv additional-output/
mv rnacode.csv additional-output/
mv rnafold-output additional-output/
mv ncrna-dfam-distance.csv additional-output/
mv conservation.csv additional-output/
mv CPC2.csv additional-output/
mv gerp.csv additional-output/
mv gnomad.csv additional-output/

#####################################################################

echo Finished calculating functionality traits.
echo

rm -rf id
rm -rf output
rm -rf test
