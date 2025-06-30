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
base_name=$(basename "${initial_data%.*}" | sed 's/-dataset//')
file_name="${output_directory}/${base_name}"

## Final Output file ##
output_file="$file_name"-intrinsic.csv      

## Temporary files ## 
<<<<<<< HEAD
output_gc="${output_directory}/${base_name}-GC.csv" 
output_complexity="${output_directory}/${base_name}-low-complexity.csv" 
=======
output_gc="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')GC.csv" 
output_dinucleotide="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dinucelotide-freq.csv" 
output_dinucleotide_tmp="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dinucelotide-tmp.csv" 
dinucleotide_seqs="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dinucelotide-seqs-tmp"
dinucleotide_seq="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')dinucelotide-seq-tmp"
>>>>>>> d1b05d496cbc2a621674d383688dd5ecad8b6c9d

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

echo "AA,AC,AG,AT,CA,CC,CpG,CT,GA,GC,GG,GT,TA,TC,TG,TT" > "$output_dinucleotide_tmp"

cat "$dinucleotide_seqs" | while read -r line
do
    echo "$line" > "$dinucleotide_seq"
    perl ./scripts/B0.1-markovProperties.pl -i "$dinucleotide_seq" -k 2 -a 'ACGT'  >> "$output_dinucleotide_tmp"
done

<<<<<<< HEAD
    # Temporary files # 
    dust_input="${file_name}-dust-input.fa"
    dust_output="${file_name}-dust-output"

    while [ "$var" -le "$last_rna_id" ]; do
    
        rna_id="RNA${var}" 

        grep -w -A 1 "$rna_id" "$initial_fasta" > "$dust_input"

        if [ -s "$dust_input" ]; then

            seq=$( grep -v "$rna_id" "$dust_input" )

            dustmasker -in "$dust_input" -out "$dust_output"            # dustmasker output: >RNAid \n start - end 

            count=$(tail -n +2 "$dust_output" | wc -l)                  # Count of low complexity regions 

            # If low sequence complexity regions recorded: calculate density in sequence; else, record zero
            if [ "$count" -gt 0 ]; then

                sum=0
                seq_len=$(echo -n "$seq" | awk '{print length}')

                tail -n +2 "$dust_output" | while read -r line; do

                    start=$(echo "$line" | awk '{print $1}')
                    end=$(echo "$line" | awk '{print $3}')
    
                    diff=$(( end - start ))

                    sum=$(( sum + diff ))

                    echo "$sum" > sum-file

                done 
        
                sum=$(cat sum-file)
                low_compl_den=$(awk "BEGIN {print $sum / $seq_len}")
        
                echo "$low_compl_den" >> "$output_complexity"

            else

                echo "0" >> "$output_complexity"
        
            fi

        else 
			echo "missing sequence in fasta file" >> "$output_complexity"
		fi

        (( var++ ))
   
    done

    rm -rf "$dust_input" "$dust_output"

fi
=======
# selected dinucleotides: GA,CpG,GG,TA

awk -F, '{print $9,$7,$11,$13}' OFS="," "$output_dinucleotide_tmp" >> "$output_dinucleotide"
>>>>>>> d1b05d496cbc2a621674d383688dd5ecad8b6c9d


## Join output files for better organization
if [ ! -s "$output_file" ]; then

    paste -d',' "$output_gc" "$output_dinucleotide" > "$output_file"

	rm -rf "$output_gc" "$output_dinucleotides"

fi


#### Remove excess files #####
<<<<<<< HEAD
#rm -rf "$output_gc" "$output_complexity"
=======
rm -rf "$dinucleotide_seqs" "$dinucleotide_seq" 
# rm -rf "$output_dinucleotide_tmp"
>>>>>>> d1b05d496cbc2a621674d383688dd5ecad8b6c9d
