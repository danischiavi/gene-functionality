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
initial_data=$1
initial_fasta=$2

first_rna_id=$(awk -F',' 'NR==2 {print $1}' "$initial_data" | tr -d 'RNA')
last_rna_id=$(awk -F',' 'END {print $1}' "$initial_data" | tr -d 'RNA') 

output_directory=data/population
mkdir -p "$output_directory"

# Final output file
output_file_variation="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')variation.csv"   

# Temporary output file
output_file_1kGP="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')1kGP-variation.csv"                  # Define name and directory for output file
output_file_gnomAD="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')gnomAD-variation.csv" 
 
# Variables
liftOver_exe=liftOver
chain_file=data/raw/hg38ToHg19.over.chain
onekGP_directory=/Volumes/archive/userdata/student_users/danielaschiavinato/dani-scratch/features-of-function-data/1000genomes
tabix_exe=tabix
vcftools_exe=vcftools 

###########################################################################################################################

# Obtain VCF Files from 1kGP 

###########################################################################################################################

#### Convert coordinates to hg19 genome version

current_coord="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')hg38-coordinates.bed"
old_coord="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')hg19-coordinates.bed"
unlifted_coord="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')unlifted.bed"

awk -F',' 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1}' "$initial_data" > $current_coord   # Reformat coordinates for LiftOver

echo "$liftOver_exe $current_coord $chain_file $old_coord $unlifted_coord &> /dev/null" >> errors.log
"$liftOver_exe" "$current_coord" "$chain_file" "$old_coord" "$unlifted_coord" &> /dev/null

##

max="$last_rna_id"
var="$first_rna_id"

## If available, unzip previously VCF files, otherwise extract files
VCF_directory_1kGP="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')1kGP-VCF" 
    
if [ -f "$VCF_directory_1kGP".zip ]; then

    unzip "$VCF_directory_1kGP".zip &> /dev/null

else
    mkdir -p "$VCF_directory_1kGP"
    
    ## Obtain VCF files according to hg19 chromosome coordinates
    while [ $var -le $max ]; do

        ## Grab the coordinates for each RNA as the counter increases
        line=$( grep -w "RNA$var" "$old_coord" | cut -f 1,2,3 | tr '\t' ' ' | tr -d "chr" | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' )
        
        ## If previous coordinate was associated with an annotated chromosome, use hg38 coordinates
        if [ -z $line ] || [[ $line == *"Un_"* ]]; then
      
            line=$( grep -w "RNA$var" $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | tr -d "chr" | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' )
            
        ## Reformatting scaffolds of annotated chromosomes
        elif [[ $line == *"_"* ]]; then
        
            a=$( echo $line | cut -d "_" -f 1 )
            b=$( echo $line | cut -d ":" -f 2 )
            line=$( echo $a:$b )

        fi
        
        rna_id='RNA'"$var"
        chr=$( echo $line | cut -d ':' -f 1 )
	
        ## The X chromosome is downloaded from a different folder to the autosomal chromosomes.
        if [ $chr == 'X' ]; then
        
            #until
	        echo "$tabix_exe -f -h "$onekGP_directory"/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz $line > $rna_id.vcf 2>>tabix.log ;"  >> tabix.log		                                   
	        #$tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
	        $tabix_exe -f -h "$onekGP_directory"/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz $line > $rna_id.vcf 2>>tabix.log ;
	        #do
	        echo Processing $rna_id.vcf >> tabix.log ; echo >> tabix.log ;
	        #./data/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	        #./data/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	        #done
            ######## Tabix is very sensitive to crashing, so repeat the command until a valid VCF file has been downloaded.

        else
            #until
	        echo "$tabix_exe -f -h "$onekGP_directory"/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz "$line" > "$rna_id".vcf 2>>tabix.log ;"  >> tabix.log
	        #$tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
	        $tabix_exe -f -h "$onekGP_directory"/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz "$line" > "$rna_id".vcf 2>>tabix.log ;
	        #do
	        echo Processing "$rna_id".vcf >> tabix.log ; echo >> tabix.log ;
	        #done
        fi                                                                            

        ######## Move VCF so it doesn't need to be redownloaded.
        mv "$rna_id".vcf "$VCF_directory_1kGP"
        
        (( var++ ))
       
        rm -rf text.file
        
    done

fi

######## Remove excess files
#rm -rf data/hg19-coordinated.bed
# rm -rf data/hg38-coordinates.bed
#rm -rf data/unlifted.bed



############################################################################################################################

# Calculate features from 1kGP

############################################################################################################################

## Function declaration
VCF2summary() {

    local input_vcf="$1"
    local len="$2"

    ## Determine SNP frequencies for each ncRNA sequence
    output_vcftools="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')vcftools-output" 
    echo "$vcftools_exe --vcf $input_vcf --freq --out $output_vcftools &> /dev/null" >> tabix.log   
    $vcftools_exe --vcf "$input_vcf" --freq --out "$output_vcftools" &> /dev/null

    ## Extract SNPs for each sequence at frequency observed at that position into a parsable format
    snps_file="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')snps"
    awk -F'\t' 'NR > 1 {print $5":"$6}' "$output_vcftools".frq > "$snps_file"    

    SNPs_count=$( wc -l < "$snps_file")
    MAF_sum=$(awk -F':' '{sum += $4} END {print sum}' "$snps_file")

    if [[ "$SNPs_count" == "0" ]] || [ -z "$SNPs_count" ]; then               # If no SNPs were recorded, set everything to NA
    
        SNPs_count='NA'
        SNPs_density='NA'
        MAF_average='NA'
    
    else

        SNPs_density=$(echo "scale=4; $SNPs_count/$len" | bc )   
	    MAF_average=$(echo "scale=4; $MAF_sum/$SNPs_count" | bc )
    
    fi
}


## Call function for each sequence
if [ ! -e "$output_file_1kGP" ]; then
    
    echo "1kGP_SNPs,1kGP_SNPSNP_density,1kGP_aveMAF" > "$output_file_1kGP"

    {
        read
        while IFS=, read -r id _ _ start end seq; do                              
    
            len=$(( $end-$start ))                                              # to calculate SNP density
            input_vcf="$VCF_directory_1kGP"/"$id".vcf
            echo "[$id][$start][$end][$len]["$id".vcf]"  >> tabix.log  

            VCF2summary "$input_vcf" "$len"                                     

            echo "$SNPs_count,$SNPs_density,$MAF_average" >> "$output_file_1kGP"
    
        done
   
    } < $initial_data   

fi

rm -rf data/hg19-coordinates.bed
rm -rf data/hg38-coordinates.bed
rm -rf data/snps
# remove VCF files or zip them? 

############################################################################################################################

# Obtain VCF Files from gnomAD and calculate features

############################################################################################################################

if [ ! -e "$output_file_gnomAD" ]; then
    
    echo "gnomAD_SNPs,gnomAD_SNP_density,gnomAD_aveMAF" > "$output_file_gnomAD"

    VCF_directory_gnomAD="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')gnomAD-VCF" 

    coordinates="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')-tabix-coordinates"

    awk -F',' 'NR > 1 {print $3":"$4"-"$5}' "$initial_data" > "$coordinates"

    ### Variables 
    max="$last_rna_id"
    var="$first_rna_id"


    while IFS= read -r line; do
        
        tabix_input="$line"

        start=$(echo "$line" | awk -F'[:-]' '{print $2}')
        end=$(echo "$line" | awk -F'[:-]' '{print $3}')                                 # For SNP density        len=$(( $end-$start ))  

        rna_id=RNA"$var"-gnomad

        chr=$(echo "$line" | awk -F'[:-]' '{print $1}')
        gnomad_database=data/raw/gnomad/gnomad.genomes.v4.0.sites."$chr".vcf.bgz               # make it a variable with user input

        echo Extracting :"$rna_id".vcf >> tabix.log ;
        echo "$tabix_exe -f -h $gnomad_database $tabix_input > $rna_id.vcf 2>>tabix.log ;" >>tabix.log
        $tabix_exe -f -h "$gnomad_database" "$tabix_input" > "$rna_id".vcf 2>>tabix.log ;
        #do
        echo >> tabix.log 
#       done

        snps_file_gnomAD="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')snps-gnomAD"
        awk -F';' '{for(i=1;i<=NF;i++) if($i ~ /^AF=/) {split($i, a, "="); print a[2]}}' "$rna_id".vcf > "$snps_file_gnomAD"
        
        SNPs_count=$( wc -l < "$snps_file_gnomAD")
        MAF_sum=$(awk -F':' '{sum += $1} END {print sum}' "$snps_file_gnomAD")

        if [[ "$SNPs_count" == "0" ]] || [ -z "$SNPs_count" ]; then               # If no SNPs were recorded, set everything to NA
    
            SNPs_count='NA'
            SNPs_density='NA'
            MAF_average='NA'
    
        else

            SNPs_density=$(echo "scale=4; $SNPs_count/$len" | bc )   
	        MAF_average=$(echo "scale=5; $MAF_sum/$SNPs_count" | bc )
    
        fi
   
        echo "$SNPs_count,$SNPs_density,$MAF_average" >> "$output_file_gnomAD" 
    
        (( var++ )) 

        mv "$rna_id".vcf "$VCF_directory_gnomAD"
    
    done < "$coordinates"

fi

if [ ! -e "$output_file_variation" ]; then

    paste -d',' "$output_file_1kGP" "$output_file_gnomAD"   > "$output_file_variation"

fi

######## Zip VCF files for further analyses
#zip -r VCF VCF &> /dev/null
#rm -rf VCF/