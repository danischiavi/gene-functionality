#!/bin/bash
#
# Script Name: B0-intrinsic-seq-features.sh
#
# Description: Calculates intrinsic sequence feature GC% and the sequence low complexity density 
#
########################################################################################################################### 

#### GENERAL SET UP ####

#### Files and directories #### 
initial_data=$1
initial_fasta=$2
             
output_directory=data/intrinsic
mkdir -p "$output_directory"
file_name="${output_directory}/$(basename "${initial_data%.*}" | sed 's/-dataset//')"

## Final Output file ##
output_file="$file_name"-intrinsic.csv      

## Temporary files ## 
output_gc="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')GC.csv" 
output_dinucleotide="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dinucelotide-freq.csv" 
output_dinucleotide_tmp="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dinucelotide-tmp.csv" 
dinucleotide_seqs="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dinucelotide-seqs-tmp"
dinucleotide_seq="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dinucelotide-seq-tmp"

########################################################################################################################### 

# GC% calculation

###########################################################################################################################

if [ ! -s "$output_gc" ]; then

    echo "GC_percentage" > "$output_gc" 
                                                       
    tail -n +2 "$initial_data"  | while IFS=, read -r _ _ _ _ _ seq; do                
    
        G_cont=$( echo "$seq" | grep -o "G\|g" | wc -l )
        C_cont=$( echo "$seq" | grep -o "C\|c" | wc -l )
        A_cont=$( echo "$seq" | grep -o "A\|a" | wc -l )
        T_cont=$( echo "$seq" | grep -o "T\|t" | wc -l )

        GC_count=$(( G_cont + C_cont ))
        total=$(( GC_count + A_cont + T_cont ))

        if (( $( echo "$GC_count == 0" | bc -l ) )); then
        
            GC=0
            percentage=NA

        else

            GC=$( echo "scale=2; $GC_count/$total" | bc )  
            percentage=$( echo "$GC*100" | bc )
        
        fi
        
        echo "$percentage" >> "$output_gc"

    done 

fi                                       

########################################################################################################################### 

# Dinucleotide frequencies 

###########################################################################################################################

# Remove header from regions file and get sequences
awk -F, 'NR > 1 {print $6}' OFS="\t", "$initial_data" > "$dinucleotide_seqs"

echo "AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT" > "$output_dinucleotide_tmp"

cat "$dinucleotide_seqs" | while read -r line
do
    echo "$line" > "$dinucleotide_seq"
    ./scripts/B0.1-markovProperties.pl -k 2 -i "$dinucleotide_seq" >> "$output_dinucleotide_tmp"
    echo "" >> "$output_dinucleotide_tmp"
done

echo "GA,CpG,GG,TA" > "$output_dinucleotide"
awk -F, '{print $9,$10,$11,$13}' OFS="," "$output_dinucleotide_tmp" > "$output_dinucleotide"


## Join output files for better organization
if [ ! -s "$output_file" ]; then

    paste -d',' "$output_gc" "$output_dinucleotide" > "$output_file"

	rm -rf "$output_gc" "$output_dinucleotides"

fi


#### Remove excess files #####
rm -rf "$dinucleotide_seqs" "$dinucleotide_seq" "$output_dinucleotide_tmp"
