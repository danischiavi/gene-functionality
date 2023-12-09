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
#        $4 is the text file containing RefSeq IDs from HGNC
#        $5 is the GFF file for the downloaded human genome Eg: GCF_000001405.39_GRCh38.p13_genomic.gff
#        $6 is the converted CSV file of the human genome. Eg: GRCh38.p13_genome.csv

###########################################################################################################################

######## General setup

###Rename parameters for better readability  #paths to files (here just for easy access) --> set them on $PATH!
#lncrna_fasta=/Volumes/archive/gardnerlab/helenacooper/features-of-function/data/lncrna/rnacentral/long-ncrna-rnacentral.fa 
#short_ncrna_fasta=/Volumes/archive/gardnerlab/helenacooper/features-of-function/data/short-ncrna/rnacentral/short-ncrna-rnacentral.fa
#chr_coord_ncrna=/Volumes/archive/gardnerlab/features_of_function/data/homo_sapiens.GRCh38.bed
# /Volumes/archive/gardnerlab/helenacooper/features-of-function/data/protein-coding/HGNC gene_with_protein_product.txt
#HGNC-protein-coding-dataset.csv

lncrna_fasta="$1"
short_ncrna_fasta="$2"
chr_coord_ncrna="$3"
ids_HGNC="$4"
genome_gff="$5"
genome_csv="$6"

# set -xv
IFS=$'\n'
d=$(date +%y%m%d) 
calc() { awk "BEGIN{print $*}"; }

#### File creation

file_creation() {

    local name=$1
    echo "ID,Functional,Chromosome,Start,End,Sequence" > ./data/$name-dataset.csv
}

names=("protein-exon2" "protein-exon3" "functional-lncrna-exon2" "functional-lncrna-exon3" "functional-short-ncrna")
for name in "${names[@]}"; do file_creation "$name"; done


######## Convert input from FASTA to CSV file for easier parsability
## check if files exist 
if [ ! -f "$lncrna_fasta" ] || [ ! -f "$short_ncrna_fasta" ] || [ ! -f "$chr_coord_ncrna" ] 
then
    echo "Initial data files do not exist"
    break 
else
    :
fi

fasta_formatter -i "$lncrna_fasta" -o "lncrna.csv" -t 
fasta_formatter -i "$short_ncrna_fasta" -o "short-ncrna.csv" -t


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

#### Arrays for searching short and lncrna ids within all ncrna RNAcentral database  

declare -a "IDs_ncrna=()"  
declare -a "IDs_lncrna=()"
declare -a "IDs_short_ncrna=()"
declare -a "IDs_pc_rna=()"

populate_array() {
    
    local file=$1
    local name=$2
    local -n array_ref="IDs_$name"

    while IFS=$'\t' read -r ID _; do
        array_ref+=("$ID")
    done < "$file"

}

populate_array lncrna.csv lncrna
populate_array short-ncrna.csv short_ncrna


while read -r _ _ _ ID _; do 
    IDs_ncrna+=("$ID")
done < "$chr_coord_ncrna"

    

#### Searching for the Id on the RNAcentral database
for id in "${IDs_ncrna[@]}"; do

    if [[ "${IDs_lncrna[@]}" =~ "$id" ]]  
    then
        set_variables "$id" "$chr_coord_ncrna" "lncrna.csv" 
        if [ "$status" != 'not-pass' ]
        then 
            exon_count=$( echo "$meta" | cut -f 10 )
            if [ "$exon_count" -ge "4" ]        #at least 4 exons
            then
                len_two=$( echo $meta | cut -f 11 | cut -d ',' -f 2 ) #length exon2 
                len_three=$( echo $meta | cut -f 11 | cut -d ',' -f 3 ) ##length exon3

                if [ "$len_two" -ge "$lower_limit" ] && [ "$len_two" -le "$upper_limit" ] || [ "$len_three" -ge "$lower_limit" ] && [ "$len_three" -le "$upper_limit" ]   # exon length cannot be greater than 3000 bp #D: ADD minimum length: -gt "75"
                then
                    ((lncrna_count++))
                    start_two=$(( $true_start + $two + 1 ))  # +1 to account for the first relative start pos being 0
                    start_three=$(( $true_start + $three + 1 )) 
                    end_two=$(( $start_two + $len_two ))
                    end_three=$(( $start_three + $len_three ))

                    echo RNA$lncrna_count,Yes,"$chr","$start_two","$end_two" >> lncrna-exon2-dataset.csv
                    echo RNA$lncrna_count,Yes,"$chr","$start_three","$end_three" >> lncrna-exon3-dataset.csv  #"ID,Functional,Chromosome,Start,End,Sequence" > $name-dataset.csv
                else
                    :
                fi
            else
                :
            fi
        fi


    elif [[ "${IDs_short_ncrna[@]}" =~ "$id" ]]
    then
        set_variables "$id" "$chr_coord_ncrna" "short-ncrna.csv" 
        if [ "$status" != 'not-pass' ]
        then
            IFS=' ' read -r chr start end <<< "$meta" #CHECK! 
            len=$(( $end - $start + 1 ))  # Add 1 to match up the sequence length to its true value????
            if [ "$len" -ge "$lower_limit" ] && [ "$len" -le "$upper_limit" ]  
            then
                ((short_count++))
                echo echo RNA$short_count,Yes,"$chr","$start","$end" >> functional-short-ncrna-dataset.csv
        
            else
                :
        
            fi 
        fi
    
    else 
        :
    fi

done




gff2Info() {

    local exons=$1
    local genome_csv=$2

    coords_one=$(   awk 'NR==1 {print $1, $4, $5}' "$exons")       # Required to generate the upstream negative control sequences
    coords_two=$(   awk 'NR==2 {print $1, $4, $5}' "$exons")       # Exon two coordinates
    coords_three=$( awk 'NR==3 {print $1, $4, $5}' "$exons")       # Exon three coordinates
    final_end=$(tail -1 "$exons" | awk '{print $1, $4, $5}')       # Required to generate the downstream negative control sequences

    chr=$( echo $coords_one | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 )   # Chromosome variable
    test=$( echo $chr | cut -c1 )   # Records any zeros in the chromosome variable
    other=$( echo $chr | cut -c2 )  # If zero is in chromosome variable, only record the single digit (ie: 01 becomes 1).
    mt_test=$( echo $coords_one | cut -d ' ' -f 1 )   # Variable to check if gene is located on the mitochondrial genome.

     # Reformat chr variable
    if [ -z "$chr" ]  # If chromosome variable empty (genes/mRNA that have been removed), then rename to allow it to be filtered out.
    then
        chr=26   
    elif [[ "$mt_test" == "NC_012920.1" ]]  # If gene is encoded on the mitochondrial genome, then rename to allow it to be filtered out.
    then
        chr=25     
    elif [[ "$test" == "0" ]] # If chromosome variable begins with zero, then rename as a single digit (ie: 01 becomes 1).
    then
        chr=$other
    elif [[ "$chr" == "23" ]]   # Chromosome X is NC_000023, but should be recorded as X in the final dataset for readability.
    then
        chr=X
    else
        :
    fi

    # Process exon data to create a dataset
    if [ "$chr" -le "22" ] || [[ "$chr" == "X" ]]   
    then
        IFS=' ' read -r _ start_one end_one     <<< "$coords_one"                                               # Starting and end coordinates of exons
        IFS=' ' read -r _ start_two end_two     <<< "$coords_two"
        IFS=' ' read -r _ start_three end_three <<< "$coords_three"
        
        seq_two=$(   grep -w "chromosome $chr" "$genome_csv" | cut -f 2 | cut -c$start_two-$end_two )           # Exons sequences 
        seq_three=$( grep -w "chromosome $chr" "$genome_csv" | cut -f 2 | cut -c$start_three-$end_three )       
        end_final=$( echo $final_end | cut -d ' ' -f 3 )                                                               # End position of final exon
	
        len_two=$(( $end_two - $start_two ))                                                                    # Length of exons
        len_three=$(( $end_three - $start_three ))          
    
        # Exclude empty and with any unknown nucleotides (N) sequences
        if [ ! -z "$seq_two" ] && [ ! -z "$seq_three" ] && [[ "$seq_two" != *"N"* ]] && [[ "$seq_three" != *"N"* ]]
        then
       
	    # Exclude sequences out of length limits 
	        if [ "$len_two" -ge "$lower_limit" ] && [ "$len_two" -le "$upper_limit" ] || [ "$len_three" -ge "$lower_limit" ] && [ "$len_three" -le "$upper_limit" ] 
	        then
               
                ((pi_count++))
                echo RNA$pi_count,Yes,chr$chr,$start_two,$end_two,$seq_two >> ./data/protein-exon2-dataset.csv
                echo RNA$pi_count,Yes,chr$chr,$start_three,$end_three,$seq_three >> ./data/protein-exon3-dataset.csv
    
                if [ "$start_one" -gt "$end_final" ]   # Reverse transcripts can alter order of start/end positions
                then
		            # File to generate negative control sequences that are the same length as exons two and three
                    echo chr$chr,$end_final,$start_one,$len_two,$len_three >> ./data/coords-for-negative-control.csv
                else
                    echo chr$chr,$start_one,$end_final,$len_two,$len_three >> ./data/coords-for-negative-control.csv
                fi
            else
                :
            fi 
        else
            :
	    fi
	else
	    :
	fi
}

## Pi coding rna

pi_count=0

while IFS= read -r id
do
    grep "exon-$id" "$genome_gff" > ./data/exons        # Grep sequence information from GFF based on RefSeq ID
    if [ "$(wc -l < ./data/exons)" -ge 4 ]; then        # At least 4 exons 
    gff2Info ./data/exons "$genome_csv"
    else
        :
    fi
done < "$ids_HGNC"