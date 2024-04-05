#!/bin/bash
#
# Script Name: intrinsic-seq-features.sh
#
# Description: This script calculates intrinsic sequence feature GC% 
#
########################################################################################################################### 

# GC% calculation

###########################################################################################################################

#### Files and directories #### 
initial_data=$1
             
output_directory=data/intrinsic
mkdir -p "$output_directory"
output_file="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')GC.csv"                  


#### Calculate GC content #### 
if [ ! -s "$output_file" ]; then

    echo "GC" > "$output_file" 
                                                       
    tail -n +2 "$initial_data"  | while IFS=, read -r _ _ _ _ _ seq; do                
    
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
        
        echo "$percentage" >> "$output_file"

    done 

fi                                       
