#!/bin/bash
#
# Script Name: slicing-database.sh
#
# Author: Daniela Schiavinato
# Last edited: November 2023
#
# Description: This script filters the functional long ncRNAs from RNAcentral database. (Future additions: filter also protein-coding and short RNA) 
#
# Input:
#        $1 is the downloaded fasta file from RNAcentral containing the lncRNAs. 
#        $2 is the downloaded chromosome coordinates file from RNAcentral for all ncRNAs in the database. Eg: Homo_sapiens.GRCh38.bed



###########################################################################################################################

######## General setup

###Rename parameters for better readability  #paths to files (here just for easy access)
lncrna_fasta=/Volumes/archive/gardnerlab/helenacooper/features-of-function/data/lncrna/rnacentral/long-ncrna-rnacentral.fa 
chr_coord_ncrna=/Volumes/archive/gardnerlab/features_of_function/data/homo_sapiens.GRCh38.bed


# set -xv
IFS=$'\n'
d=$(date +%y%m%d) 

#######Declare function

function my_function() {

    local id=$1
    local chr_coord=$2
    local seq_csv=$3
    local n=$4
    
    meta=$( grep -m 1 "$id" "$chr_coord" )
    chr=$( echo "$meta" | cut -f 1 )
    true_start=$( echo "$meta" | cut -f 2 )  
    true_end=$( echo "$meta" | cut -f 3 )
    len=$( echo "$meta" | cut -f 11 | cut -d ',' -f 1 )
    end=$(( $true_start + $len ))
    exon_count=$( echo "$meta" | cut -f 10 )
    seq=$( grep -m 1 "$id" "$seq_csv" | cut -f 2 )
        
    if [ -z "$seq" ]   # If no sequence available for lncRNA, remove
    then
        :

	elif [ $chr == 'chrM' ] || [ $chr == 'chrY' ]  # Removes ncRNA from mitochondria and Y chromosome
    then
        :
	
    elif [ $exon_count -lt "4" ]   # Need at least four exons
    then
        :
	    
    else
        len_two=$( echo $meta | cut -f 11 | cut -d ',' -f 2 ) #length exon2 
        len_three=$( echo $meta | cut -f 11 | cut -d ',' -f 3 ) ##length exon3

        if [ "$len_two" -lt "3000" ] || [ "$len_three" -lt "3000" ]   # exon length cannot be greater than 3000 bp #D: ADD minimum length: -gt "75"
        then
            echo $meta >> $d-$n-sliced-database.bed
        else
            :
        fi
    fi
}

 
######## Convert input from FASTA to CSV file for easier parsability
fasta_formatter -i "$lncrna_fasta" -o "$d-converted-lncrna.csv" -t 


######## Obtain chromosome coordinates for each functional ncRNA

####Define assosiative array for more efficient searching

IDs_ncrna=()  #to store the IDs of all ncRNA from RNAcentral database 
IDs_long_ncrna=()   #to store the IDs of long-ncRNA from RNAcentral database


##Populate the arrays

while read -r line; do
    ID=$( echo "$line" | cut -f 1 | cut -c1-18 )
    IDs_long_ncrna+=("$ID")
done < "$d-converted-lncrna.csv"

while read -r line; do 
    IdField=$( echo "$line" | cut -f 4 | cut -c1-18 )
    IDs_ncrna+=("$IdField")
done < "$chr_coord_ncrna"


####Searching for the Id on the RNAcentral database
for id in "${IDs_ncrna[@]}"; do

    id_var="$id"

    if [[ "${IDs_long_ncrna[@]}" =~ "$id_var" ]]
    then
        my_function "$id_var" "$chr_coord_ncrna" "$d-converted-lncrna.csv" "lncrna"

    else
        :
    fi 

done