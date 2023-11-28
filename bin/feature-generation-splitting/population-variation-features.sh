#!/bin/bash
#
# Script Name: population-variation-features.sh
#
#
# Description: This script calculates the chosen population variation features for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences including SNP count, SNP density and minor allele frequency. 
#
# Input: $1 is the the dataset
#        $2 is the location of the folder with the local databases/datasets or version specific executables
# 
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.

############################################################################################################################

#### Variable declarations
initial_data=$1
additional_folder=$2 

if find $additional_folder/ -executable -type f | grep -q vcftools #maybe this could be in the runall script? 
then
    vcftools_exe=$additional_folder/vcftools
elif command -v vcftools &> /dev/null
then
    vcftools_exe=vcftools
else
    echo Please check that vcftools is installed.
    #exit 1
fi


#### File creation 
echo 1000G_SNPs,1000G_SNPsDensity,aveMAF > 1000g-freqsummary.csv


##### Function declarations
VCF2summary() {

    local vfctools_exe="$1"
    local vcf_input="$2"
    local output_prefix="$3"
    local len="$4"

    ######## Determine SNP frequencies for each ncRNA sequence
    #echo "$vcftools_exe --vcf $input_vcf --freq --out $output_prefix &> /dev/null" >> tabix.log  #are they doing the same? 
    $vcftools_exe --vcf $input_vcf --freq --out $output_prefix &> /dev/null

    ######## Extract each SNP at frequency observed at that position into a parsable format
    grep -v "CHROM" $output_prefix.frq | cut -f 5,6 | perl -lane '{print "$F[0]:$F[1]"}' > snps.file  
    min=0  # Min MAF observed?? it was min=0.99, but I think should be min=0! 
    sum=0
    count=0

    while IFS=: read -r snpa _ snpb maf  # # Add a way to Stop loop crashing if VCF file is empty
    do
        sum=$(echo "$sum + $maf" | bc )

        # I think all this section is not doing anything?? 
        ######## If the maximum? (dont really understand what max is.. I think it should be minimun?? ) is recorded as zero, set the new values automatically as the max D:???? where is max defined? 
        #if (( $(echo "$min == 0" | bc -l) ))
        #then
        #   min=$maf
        ######## If the new SNP analysed has a frequency smaller than the recorded min, update the min
        #elif (( $(echo "$maf < $min" | bc -l) ))
        #then
        #    min=$maf
        #else
        #    :
        #fi
        

        (( count++ ))

    done < snps.file

    if [[ "$count" == "0" ]] || [ -z "$count" ]   # If no SNPs were recorded, set everything to NA
    then
        count='NA'
        density='NA'
        average='NA'
    else
        density=$(echo "scale=3; $count/$len" | bc )   # SNP Density
	    average=$(echo "scale=3; $sum/$count" | bc )
    fi


}


######## Loop through VCF files to calculate features

{
    read
    while IFS=, read -r id _ _ start end seq
    do
        len=$(( $end-$start ))   # Sequence length
        f="$id-1k.vcf"
        input_vcf="VCF/$f"
        output_prefix="freq-afr"
        echo "[$id][$start][$end][$len][$f]"  >> tabix.log  

        VCF2summary "$vcftools_exe" "$input_vcf" "$output_prefix" "$len"

        echo "$count,$density,$average" >> 1000g-freqsummary.csv
    
    done
   
} < $initial_data


