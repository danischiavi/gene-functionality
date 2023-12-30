#!/bin/bash
#
# Script Name: intrinsic-seq-features.sh
#
# Description: This script calculates intrinsic sequence feature GC% 
#
########################################################################################################################### 

# GC% calculation

###########################################################################################################################

echo GC_percentage > data/GC-dataset-$name.csv            # Global variable defined in runall.sh for dataset name  

{
    read                                                    #skips .csv header 
    while IFS=, read -r _ _ _ _ _ seq                       #6th field .csv: sequence
    do
        G_cont=$( echo "$seq" | grep -o "G\|g" | wc -l )
        C_cont=$( echo "$seq" | grep -o "C\|c" | wc -l )
        A_cont=$( echo "$seq" | grep -o "A\|a" | wc -l )
        T_cont=$( echo "$seq" | grep -o "T\|t" | wc -l )

        GC_count=$(( $G_cont + $C_cont ))
        total=$(( $GC_count + $A_cont + $T_cont ))

        if (( $( echo "$GC_count == 0" | bc -l ) ))
        then
            GC=0
            percentage=NA
        else
            GC=$( echo "scale=2; $GC_count/$total" | bc )  #D:Check standard digits?
            percentage=$( echo "$GC*100" | bc )
        fi

        echo $percentage >> data/GC-dataset-$name.csv 
    done

} < "$initial_data"                                        # Global variable defined in runall.sh for dataset reference
