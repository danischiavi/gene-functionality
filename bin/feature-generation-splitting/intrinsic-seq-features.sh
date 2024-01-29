#!/bin/bash
#
# Script Name: intrinsic-seq-features.sh
#
# Description: This script calculates intrinsic sequence feature GC% 
#
########################################################################################################################### 

# GC% calculation

###########################################################################################################################

initial_data=$1
output_file="${initial_data/dataset.csv/GC.csv}"            # Define name of output file
echo GC_percentage > $output_file              

{
    read                                                    # Skips .csv header 
    while IFS=, read -r _ _ _ _ _ seq; do                   # 6th field .csv: sequence
    
        G_cont=$( echo "$seq" | grep -o "G\|g" | wc -l )
        C_cont=$( echo "$seq" | grep -o "C\|c" | wc -l )
        A_cont=$( echo "$seq" | grep -o "A\|a" | wc -l )
        T_cont=$( echo "$seq" | grep -o "T\|t" | wc -l )

        GC_count=$(( $G_cont + $C_cont ))
        total=$(( $GC_count + $A_cont + $T_cont ))

        if (( $( echo "$GC_count == 0" | bc -l ) )); then
        
            GC=0
            percentage=NA

        else

            GC=$( echo "scale=2; $GC_count/$total" | bc )  
            percentage=$( echo "$GC*100" | bc )
        fi
        
        echo $percentage >> $output_file

    done

} < "$initial_data"                                         
