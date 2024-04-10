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
output_complexity="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')low-complexity.csv" 

########################################################################################################################### 

# GC% calculation

###########################################################################################################################

if [ ! -s "$output_gc" ]; then

    echo "GC_percent" > "$output_file" 
                                                       
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

# Low complexity density 

###########################################################################################################################

if [ ! -s "$output_complexity" ]; then

    echo "lowComplexity_density" > "$output_complexity"

    # Variables to read fasta file # 
    first_rna_id=$(awk 'NR==1 {print $1}' "$initial_fasta" | tr -d '>RNA') 
    last_rna_id=$(awk '/^>/{id=$1} END{print id}' "$initial_fasta" | tr -d '>RNA')
    var="$first_rna_id"

    # Temporary files # 
    dust_input="$file_name"-dust-input.fa
    dust_output="$file_name"-dust-output

    while [ "$var" -le "$last_rna_id" ]; do
    
        rna_id="RNA$var" 

        grep -w -A 1 "$rna_id" "$initial_fasta" > "$dust_input"

        if [ -s "$dust_input" ]; then

            seq=$( grep -v "$rna_id" "$dust_input" )

            dustmasker -in "$dust_input" -out "$dust_output"            # dustmasker output: >RNAid \n start - end 

            count=$(tail -n +2 "$dust_output" | wc -l)                  # Count of low complexity regions 

            # If low sequence complexity regions recorded: calculate density in sequence; else, record zero
            if [ "$count" -gt 0 ]; then

                sum=0
                seq_len=$(echo -n "$seq" | awk '{print length}')

                tail -n +2 TEST1 | while read -r line; do

                    start=$(echo "$line" | awk '{print $1}')
                    end=$(echo "$line" | awk '{print $3}')
    
                    diff=$((end - start))

                    sum=$((sum + diff))

                    echo "$sum" > sum-file

                done 
        
                sum=$(cat sum-file)
                low_compl_den=$(awk "BEGIN {print $sum / $seq_len}")
        
                echo "$low_compl_den" >> "$output_complexity"

            else

                echo "0" >> "$output_complexity"
        
            fi
        fi

        (( var++ ))
   
    done
fi


## Join output files for better organization
if [ ! -s "$output_file" ]; then

    paste -d',' "$output_gc" "$output_complexity" > "$output_file"

fi


#### Remove excess files #####
#rm -rf "$output_gc"
#rm -rf "$output_complexity"