#!/bin/bash
#
# Script Name: population-variation-features.sh
#
#
# Description: This script calculates the chosen population variation features for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences including SNP count, SNP density and minor allele frequency. 
#
# Input: $1 is the the dataset
#        $2 is the location of the folder with the local gnomAD databases
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.

############################################################################################################################
##General Set up ## 
initial_data=$1
gnomad_directory=$2

first_rna_id=$(awk -F',' 'NR==2 {print $1}' "$initial_data" | tr -d 'RNA') 

output_directory=data/population
mkdir -p "$output_directory"

# Final output file
file_name=${output_directory}/$(basename "${initial_data%.*}" | sed 's/-dataset//')
output_file="$file_name"-population.csv 

# Variables
tabix_exe=tabix

############################################################################################################################

# Obtain VCF Files from gnomAD and calculate features

############################################################################################################################
VCF_directory_gnomAD="$file_name"-gnomAD-VCF

if [ ! -s "$output_file" ]; then
    
    echo "SNP_density,MAF_avg" > "$output_file"
    
    if [ ! -e "$VCF_directory_gnomAD" ]; then
    
        if [ -f "$VCF_directory_gnomAD".zip ]; then

            unzip "$VCF_directory_gnomAD".zip &> /dev/null
        
        else

            mkdir -p "$VCF_directory_gnomAD"

        fi

    fi

    ### Variables ##
    coordinates="$file_name"-tabix-coordinates
    awk -F',' 'NR > 1 {print $3":"$4"-"$5}' "$initial_data" > "$coordinates"

    var="$first_rna_id"

    ## Calculate SNPs ## 
    while IFS= read -r line; do
        
        rna_id="$VCF_directory_gnomAD"/RNA"$var"
        
        if [ ! -e "$rna_id" ]; then

            tabix_input="$line"

            start=$(echo "$line" | awk -F'[:-]' '{print $2}')
            end=$(echo "$line" | awk -F'[:-]' '{print $3}')                                         
            len=$(( $end - $start + 1 ))                                                               # For SNP density. +1 since coordinates are inclusive  

            chr=$(echo "$line" | awk -F'[:-]' '{print $1}')
            gnomad_database="$gnomad_directory"/gnomad.genomes.v4.0.sites."$chr".vcf.bgz               # make it a variable with user input

            echo "Extracting:$rna_id.vcf" >> tabix.log 
            echo "$tabix_exe -f -h $gnomad_database $tabix_input > $rna_id.vcf 2>>tabix.log" >>tabix.log
            
            $tabix_exe -f -h "$gnomad_database" "$tabix_input" > "$rna_id".vcf 2>>tabix.log 
            #do
            echo >> tabix.log 
#           done
        fi

        snps_file="$file_name"-snps
        awk -F';' '{for(i=1;i<=NF;i++) if($i ~ /^AF=/) {split($i, a, "="); print a[2]}}' "$rna_id".vcf > "$snps_file"
        
        SNPs_count=$( wc -l < "$snps_file")
        MAF_sum=$(awk -F':' '{sum += $1} END {printf "%.10f", sum}' "$snps_file")
      

        if [[ "$SNPs_count" == "0" ]] || [ -z "$SNPs_count" ]; then                                     # If no SNPs were recorded, set everything to NA
    
            SNPs_density='NA'
            MAF_average='NA'
    
        else

            SNPs_density=$(echo "scale=6; $SNPs_count/$len" | bc )   
	        MAF_average=$(echo "scale=6; $MAF_sum/$SNPs_count" | bc )
    
        fi

        echo "$SNPs_density,$MAF_average" >> "$output_file" 
    
        (( var++ )) 
    
    done < "$coordinates"

    ## Zip VCF files for further analyses
    zip -r "$VCF_directory_gnomAD".zip "$VCF_directory_gnomAD" &> /dev/null && rm -rf "$VCF_directory_gnomAD" 

fi

### Remove excess files
rm -rf "$snps_file" 
rm -rf "$coordinates"
