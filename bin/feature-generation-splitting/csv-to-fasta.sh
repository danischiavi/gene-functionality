#!/bin/bash

csv_to_fasta() {                     
    
    local file=$1
    {
    read 
    while IFS=',' read -r id _ _ _ _ seq; do
            output_file="${file/dataset.csv/seq.fa}"
            echo -e ">$id\n$seq" >> $output_file
        
    done

    } < $file   
}   

csv_to_fasta $1