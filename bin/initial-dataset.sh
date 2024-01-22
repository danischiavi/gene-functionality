#!/bin/bash
#
# This script is equivalent to slicing-database but implementing Snakemake --> see Snakefile 


###########################################################################################################################

# General setup

###########################################################################################################################

# Inputs
rnacentral_coords= "{input.rnacentral_coords}"
rnacentral_lncrna="{input.rnacentral_lncrna}"
rnacentral_short_ncrna="{input.rnacentral_short_ncrna}"
rnacentral_pre_mirna="{input.rnacentral_pre_mirna}"
genome_annotations="{input.genome_annotations}"
genome_seq="{input.genome_seq}"
protein_coding_refseq="{input.protein_coding_refseq}"

# Outputs
lncrna_exon_two_csv="{output.lncrna_exon_two_csv}"
lncrna_exon_three_csv="{output.lncrna_exon_three_csv}"
negative_control_csv="{output.negative_control_csv}"
short_ncrna_csv="{output.short_ncrna_csv}"
pre_mirna_csv="{output.pre_mirna_csv}"
protein_exon_two_csv="{output.protein_exon_two_csv}"
protein_exon_three_csv="{output.protein_exon_three_csv}"


## check if files exist 
if [ ! -f "$lncrna_fasta" ] || [ ! -f "$short_ncrna_fasta" ] || [ ! -f "$chr_coord_ncrna" ]  ## ADD THE OTHER FILES!
then
    echo "Initial data files do not exist"
    break 
else
    :
fi

## File creation
file_creation() {

    local name=$1
    echo "ID,Functional,Chromosome,Start,End,Sequence" > data/$name-dataset.csv
}

names=("protein-exon2" "protein-exon3" "functional-lncrna-exon2" "functional-lncrna-exon3" "functional-short-ncrna" "functional-pre-mirna")
for name in "${names[@]}"; do file_creation "$name"; done


## Sequence length limits
lower_limit='75'
upper_limit='3000'


###########################################################################################################################

# NON-CODING RNA - Chromosome coordinates and sequence 

###########################################################################################################################

#### Function declaration

populate_array() {                                ## Populate arrays with short and lncrna ids for searching within all ncrna RNAcentral database
    
    local file=$1
    local name=$2
    local -n array_ref="IDs_$name"

    while IFS=$'\t' read -r ID_full _; do         # ID with description (rnacentral) 
       IFS=' ' read -r ID _     <<< "$ID_full"
       array_ref+=("$ID")
    done < "$file"
}


set_variables() {                                  ## Set variables to filter out ncRNA sequences from the RNAcentral database
    
    local id=$1
    local seq_csv=$2
    
    seq=$(  grep -m 1 "$id" "$seq_csv" | cut -f 2 )
    meta=$( grep -m 1 "$id" "$rnacentral_coords" )
    chr=$(  echo "$meta" | cut -f 1 )
    status='pass'                                  # Set status to filter out sequences which are redifine as 'not-pass'
        
    if [ -z "$seq" ]                               # Remove if no sequence available
    then
        status='not-pass'
        return
	elif [ $chr == 'chrM' ] || [ $chr == 'chrY' ]  # Remove ncRNA from mitochondria and Y chromosome
    then
        status='not-pass'
        return
    else
        :
    
    fi
}
                

## Arrays for searching short and lncrna ids within all ncrna RNAcentral database  

declare -a "IDs_ncrna=()"  
declare -a "IDs_lncrna=()"
declare -a "IDs_short_ncrna=()"
declare -a "IDs_pre_mirna=()"

# Populate lncRNA and short-ncRNA arrays
populate_array $rnacentral_lncrna "lncrna"
populate_array $rnacentral_short_ncrna "short_ncrna"
populate_array $rnacentral_pre_mirna "pre_mirna"

# Populate ncRNA array
while IFS=$'\t' read -r _ _ _ ID _; do 
    if [[ "$ID" != "${IDs_ncrna[@]}" ]]; then
        IDs_ncrna+=("$ID")
    else  
        :
    fi
done < "$rnacentral_coords"
   
## D: ADD CHECKS FOR ARRAYS!    

## Searching for the Id on the RNAcentral database
lncrna_count=0
short_count=0
pre_mirna_count=0

for id in "${IDs_ncrna[@]}"; do

    if [[ "${IDs_lncrna[@]}" =~ "$id" ]]  
    then
        set_variables "$id" $rnacentral_lncrna 
        if [ "$status" != 'not-pass' ]                                                                      # Status returned by function if empty sequence or mitocondrial/Ychr
        then 
            exon_count=$( echo "$meta" | cut -f 10 )
            if [ "$exon_count" -ge 3 ]                                                                      # Filter for sequences with at least 4 exons
            then 
                len_two=$(  echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $2}')               # Length exons within range
                len_three=$(echo "$meta" | awk -F'\t' '{print $11}' | awk -F',' '{print $3}') 
                len_last=$( echo "$meta" | awk -F'\t' '{print $11}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')

                if ([ "$len_two" -ge "$lower_limit" ] && [ "$len_two" -le "$upper_limit" ]) && ([ "$len_three" -ge "$lower_limit" ] && [ "$len_three" -le "$upper_limit" ])  
                then
                    ((lncrna_count++))

                    seq_start=$(echo "$meta" | awk -F'\t' '{print $2}' )                                     # 0-start. Function in progress for this (see ./bin/draft)
                    relative_start_one=$(  echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $1}') # Position relative to Start (bed format) for exons     
                    relative_start_two=$(  echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $2}')     
                    relative_start_three=$(echo "$meta" | awk -F'\t' '{print $12}' | awk -F',' '{print $3}')
                    relative_start_last=$( echo "$meta" | awk -F'\t' '{print $12}' | awk -v exon_count="$exon_count" -F',' '{print $exon_count}')  
                    
                    start_one=$((   $seq_start + $relative_start_one + 1 ))                                  # +1 to account for the first relative start pos being 0
                    start_two=$((   $seq_start + $relative_start_two + 1 ))                                         
                    start_three=$((              $seq_start + $three + 1 ))
                    start_last=$(( $seq_start + $relative_start_last + 1 )) 
                    end_two=$((                    $start_two + $len_two ))
                    end_three=$((              $start_three + $len_three ))
                    end_last=$((                 $start_last + $len_last ))
                                    
                    echo RNA$lncrna_count,Yes,"$chr","$start_two","$end_two","$seq"     >> "$lncrna_exon_two_csv"
                    echo RNA$lncrna_count,Yes,"$chr","$start_three","$end_three","$seq" >> "$lncrna_exon_three_csv"  

                    if [ "$start_one" -gt "$end_last" ]                                                       # Reverse transcripts can alter order of start/end positions
                    then
		            # To generate negative control sequences that are the same length as exons two and three
                        echo $chr,$end_last,$start_one,$len_two,$len_three >> "$negative_control_csv"   
                    else
                        echo $chr,$start_one,$end_last,$len_two,$len_three >> "$negative_control_csv"
                    fi
                fi
            fi
        fi


    elif [[ "${IDs_short_ncrna[@]}" =~ "$id" ]]
    then
        set_variables "$id" "$rnacentral_short_ncrna" 
        if [ "$status" != 'not-pass' ]
        then
            IFS=$'\t ' read -r chr zero_start end _ <<< "$meta"                                              # zero_start: 0-start bed format
            len=$(( $end-$zero_start ))  
            if [ "$len" -ge "$lower_limit" ] && [ "$len" -le "$upper_limit" ]  
            then
                (( short_count++ ))
                start=$((zero_start+=1))                                                                     # add 1 to change coordinate from 0-start
                echo RNA$short_count,Yes,"$chr","$start","$end","$seq" >> "$short_ncrna_csv"
                
                if [ "$start" -gt "$end" ]                                                                    # Reverse transcripts can alter order of start/end positions
                then
		            # To generate negative control sequences that are the same length as exons two and three
                    echo $chr,$end,$start,$len,'NA' >> "$negative_control_csv"
                else
                    echo $chr,$start,$end,$len,'NA' >> "$negative_control_csv"
                fi
            fi 
        fi
    elif [[ "${IDs_pre_mirna[@]}" =~ "$id" ]]
    then
        set_variables "$id" "$rnacentral_pre_mirna" 
        if [ "$status" != 'not-pass' ]
        then
            IFS=$'\t ' read -r chr zero_start end _ <<< "$meta"  
            len=$(( $end - $zero_start ))  
            if [ "$len" -ge "$lower_limit" ] && [ "$len" -le "$upper_limit" ]  
            then
                (( pre_mirna_count++ ))
                start=$((zero_start+=1)) 
                echo RNA$pre_mirna_count,Yes,"$chr","$start","$end","$seq" >> "$pre_mirna_csv"
                
                if [ "$start" -gt "$end" ]                                                                     # Reverse transcripts can alter order of start/end positions
                then
		            # To generate negative control sequences that are the same length as exons two and three
                    echo $chr,$end,$start,$len,'NA' >> "$negative_control_csv"
                else
                    echo $chr,$start,$end,$len,'NA' >> "$negative_control_csv"
                fi
            fi 
        fi   
    fi
done


###########################################################################################################################

# PROTEIN CODING RNA - Chromosome coordinates and sequence 

###########################################################################################################################

##### Function declaration ####

gff2Info() {                                                                                ## Get info from the Human genome database gff file to filter out pi-coding sequences 

    local exons=$1
    local genome_csv=$2

    coords_one=$(   awk 'NR==1 {print $1, $4, $5}' "$exons")                                # Required to generate the upstream negative control sequences
    coords_two=$(   awk 'NR==2 {print $1, $4, $5}' "$exons")                                # Exon two coordinates
    coords_three=$( awk 'NR==3 {print $1, $4, $5}' "$exons")                                # Exon three coordinates
    final_end=$(tail -1 "$exons" | awk '{print $1, $4, $5}')                                # Required to generate the downstream negative control sequences

    chr=$( echo $coords_one | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 ) # Chromosome variable
    test=$(  echo $chr | cut -c1 )                                                          # Records any zeros in the chromosome variable
    other=$( echo $chr | cut -c2 )                                                          # If zero is in chromosome variable, only record the single digit (ie: 01 becomes 1).
    mt_test=$( echo $coords_one | cut -d ' ' -f 1 )                                         # Variable to check if gene is located on the mitochondrial genome.

    # Reformat chr variable or rename to allow it to be filtered out
    if [ -z "$chr" ]                                                                        # If chromosome variable empty (genes/mRNA that have been removed)
    then
        chr=26   
    elif [[ "$mt_test" == "NC_012920.1" ]]                                                  # If gene is encoded on the mitochondrial genome
    then
        chr=25     
    elif [[ "$test" == "0" ]]                                                               # If chromosome variable begins with zero, then rename as a single digit (ie: 01 becomes 1).
    then
        chr="$other"    
    fi

    # Process exon data to create a dataset
    if [ "$chr" -le '23' ]  
    then
        IFS=' ' read -r _ start_one end_one     <<< "$coords_one"                                               # Starting and end coordinates of exons
        IFS=' ' read -r _ start_two end_two     <<< "$coords_two"
        IFS=' ' read -r _ start_three end_three <<< "$coords_three"
        
        if [[ "$chr" == 23 ]]; then chr=X; else :; fi                                                           # Chromosome X is NC_000023, but should be recorded as X in the final dataset for readability.
        
        seq_two=$(   grep -w "chromosome $chr" "$genome_csv" | cut -f 2 | cut -c$start_two-$end_two )           # Exons sequences 
        seq_three=$( grep -w "chromosome $chr" "$genome_csv" | cut -f 2 | cut -c$start_three-$end_three )       
        end_final=$( echo $final_end | cut -d ' ' -f 3 )                                                        # End position of final exon
	
        len_two=$(( $end_two - $start_two ))                                                                    # Length of exons
        len_three=$(( $end_three - $start_three ))          
        
        # Exclude empty and with any unknown nucleotides (N) sequences
        if [ ! -z "$seq_two" ] && [ ! -z "$seq_three" ] && [[ "$seq_two" != *"N"* ]] && [[ "$seq_three" != *"N"* ]]
        then
       
	    # Exclude sequences out of length limits 
	        if ([ "$len_two" -ge "$lower_limit" ] && [ "$len_two" -le "$upper_limit" ]) && ([ "$len_three" -ge "$lower_limit" ] && [ "$len_three" -le "$upper_limit" ]) 
	        then
                ((pi_count++))
                echo RNA$pi_count,Yes,chr$chr,$start_two,$end_two,$seq_two >> $protein_exon_two_csv
                echo RNA$pi_count,Yes,chr$chr,$start_three,$end_three,$seq_three >> $protein_exon_three_csv
    
                if [ "$start_one" -gt "$end_final" ]                                                            # Reverse transcripts can alter order of start/end positions
                then
		            # To generate negative control sequences that are the same length as exons two and three
                    echo chr$chr,$end_final,$start_one,$len_two,$len_three >> $negative_control_csv
                else
                    echo chr$chr,$start_one,$end_final,$len_two,$len_three >> $negative_control_csv
                fi
            fi 
	    fi
	fi
}

###############################

# Apply to each sequence

protein_count=0

while IFS= read -r id
do
    grep "exon-$id" "$genome_annotations" > data/exons        # Grep annotation from Reference Genome (NCBI) according to protein-coding genes (HGNC)
    if [ "$(wc -l < data/exons)" -ge 4 ]; then                # At least 4 exons 
    gff2Info data/exons "$genome_seq"
    fi
done < "$protein_coding_refseq"



