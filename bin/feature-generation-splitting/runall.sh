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

#### Variable declarations
initial_data=$1-dataset.csv 
initial_fasta=$1-seq.fa 

if [ ! -f $initial_data ] || [ ! -f $initial_fasta ]
then
    echo "Initial data files do not exist."
    exit 1
else
    :
fi

additional_folder=$2 
if [ ! -d $additional_folder/ ]
then
    echo "Folder with required local databases does not exist."
    #exit 1
else
    :
fi

#### Reformat initial dataset
awk -F',' 'NR > 1 {print $3, $4, $5}' $initial_data > coordinates
awk -F',' 'NR > 1 {print $3":"$4"-"$5}' $initial_data > locations
awk -F',' 'NR > 1 {print $3"\t"$4"\t"$5}' "$initial_data" | sort -k1,1 -k2,2n > sorted-coordinates.bed
awk -F',' 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1}' "$initial_data" > coordinates-rnaid.bed


############################################################################################################################

# Calculate features

############################################################################################################################

#do I have to add the files as variables for the scripts? 

./intrinsic-seq-features.sh $initial_data

./seq-conservation-features.sh coordinates $additional_folder #with $ or without? 

./transcriptome-expression-features.sh locations $additional_folder

./genomic-repeat-associated-features.sh $initial_data $initial_fasta sorted-coordinates.bed $additional_folder

./protein-and-rna-specific-features.sh xxx

./population-variation-features.sh coordinates-rnaid.bed $initial_data $additional_folder 


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
