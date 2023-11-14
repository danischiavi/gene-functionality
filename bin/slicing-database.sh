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
#        $2 is the downloaded fasta file from RNAcentral containing the short-ncRNAs.
#        $3 is the downloaded chromosome coordinates file from RNAcentral for all ncRNAs in the database. Eg: Homo_sapiens.GRCh38.bed


###########################################################################################################################

######## General setup

###Rename parameters for better readability  #paths to files (here just for easy access) --> set them on $PATH!
#lncrna_fasta=/Volumes/archive/gardnerlab/helenacooper/features-of-function/data/lncrna/rnacentral/long-ncrna-rnacentral.fa 
#short_ncrna_fasta=/Volumes/archive/gardnerlab/helenacooper/features-of-function/data/lncrna/rnacentral/short-ncrna-rnacentral.fa 
#chr_coord_ncrna=/Volumes/archive/gardnerlab/features_of_function/data/homo_sapiens.GRCh38.bed


lncrna_fasta="$1"
short_ncrna_fasta="$2"
chr_coord_ncrna="$3"

######## Convert input from FASTA to CSV file for easier parsability
## check if files exist 
if [ ! -f "$lncrna_fasta" ] || [ ! -f "$short_ncrna_fasta" ] || [ ! -f "$chr_coord_ncrna" ] 
then
    echo "Initial data files do not exist"
    break 
else
    :
fi

fasta_formatter -i "$lncrna_fasta" -o "$d-converted-lncrna.csv" -t 
fasta_formatter -i "$short_ncrna_fasta" -o "$d-converted-short-ncrna.csv" -t

# set -xv
IFS=$'\n'
d=$(date +%y%m%d) 


## sequence length limits
lower_limit='75'
upper_limit='3000'

####### Declare function

function set_variables() { 
    local id=$1
    local chr_coord=$2
    local seq_csv=$3
    
    seq=$( grep -m 1 "$id" "$seq_csv" | cut -f 2 )
    meta=$( grep -m 1 "$id" "$chr_coord" )
    chr=$( echo "$meta" | cut -f 1 )
    status='pass'
        
    if [ -z "$seq" ]   # If no sequence available for lncRNA, remove
    then
        status='not-pass'
        return
	elif [ $chr == 'chrM' ] || [ $chr == 'chrY' ]  # Removes ncRNA from mitochondria and Y chromosome
    then
        status='not-pass'
        return
    else
        :
    
    fi
}

 
######## Obtain chromosome coordinates for each functional ncRNA

####Define assosiative array for more efficient searching

IDs_ncrna=()  #to store the IDs of all ncRNA from RNAcentral database 
IDs_lncrna=()   #to store the IDs of long-ncRNA from RNAcentral database
IDs_short_ncrna=()

##Populate arrays

while read -r line; do
    ID=$( echo "$line" | cut -f 1 | cut -c1-18 )  #D: check how long the field is 
    IDs_lncrna+=("$ID")
done < "$d-converted-lncrna.csv"

while read -r line; do
    ID=$( echo "$line" | cut -f 1 | cut -c1-18 ) 
    IDs_short_ncrna+=("$ID")
done < "$d-converted-short.csv"

while read -r line; do 
    IdField=$( echo "$line" | cut -f 4 | cut -c1-18 )
    IDs_ncrna+=("$IdField")
done < "$chr_coord_ncrna"


####Searching for the Id on the RNAcentral database
for id in "${IDs_ncrna[@]}"; do

    if [[ "${IDs_lncrna[@]}" =~ "$id" ]]  #D: change order? 
    then
        set_variables "$id" "$chr_coord_ncrna" "$d-converted-lncrna.csv" 
        if [ "$status" != 'not-pass' ]
        then 
            exon_count=$( echo "$meta" | cut -f 10 )
            if [ "$exon_count" -lt "4" ]
            then
                :
        
            else 
                len_two=$( echo $meta | cut -f 11 | cut -d ',' -f 2 ) #length exon2 
                len_three=$( echo $meta | cut -f 11 | cut -d ',' -f 3 ) ##length exon3

                if [ "$len_two" -ge "$lower_limit" ] && [ "$len_two" -le "$upper_limit" ] || [ "$len_three" -ge "$lower_limit" ] && [ "$len_three" -le "$upper_limit" ]   # exon length cannot be greater than 3000 bp #D: ADD minimum length: -gt "75"
                then
                    echo "$meta" >> $d-lncrna-sliced-database-TEST.bed
                else
                    :
                fi
            fi
        fi


    elif [[ "${IDs_short_ncrna[@]}" =~ "$id" ]]
    then
        set_variables "$id" "$chr_coord_ncrna" "$d-converted-short-ncrna.csv" 
        if [ "$status" != 'not-pass' ]
        then
            len=$( echo "$meta" | cut -f 11 | cut -d ',' -f 1 )
            if [ "$len" -ge "$lower_limit" ] && [ "$len" -le "$upper_limit" ]  
            then
                echo "$meta" >> $d-short-ncrna-sliced-database-TEST.bed
        
            else
                :
        
            fi 
        fi
    fi

done







