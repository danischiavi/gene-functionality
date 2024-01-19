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

####Checks 

##vcftools
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

## tabix
if find $additional_folder/ -executable -type f | grep -q tabix
then
    tabix_exe=$additional_folder/tabix
elif command -v tabix &> /dev/null
then
    tabix_exe=tabix
else
    echo Please check that tabix is installed.
    exit 1
fi

## gnomAD VCF file

if [ -f $additional_folder/gnomad/gnomad.genomes.r3.0.sites.vcf.bgz ] 
then
    gnomad_database=$additional_folder/gnomad/gnomad.genomes.r3.0.sites.vcf.bgz
else
    until [ -f $gnomad_database ] ; do read -p "Please enter custom name of gnomAD SNP database (vcf.bz): " i ; gnomad_database=$addition
al_folder/$i ; echo ; done
fi

######## Conversion file for changing coordinates between hg19 and hg38 for liftOver

if [ -f $additional_folder/hg38ToHg19.over.chain ]
then
    chain_file=$additional_folder/hg38ToHg19.over.chain
else
    until [ -f $chain_file ] ; do read -p "Please enter custom name of coordinates file for converting between hg38 and Hg19 genome versions (chain): " chain ; chain_file=$additional_file/$chain ; echo ; done
fi

######## liftOver

if find $additional_folder/ -executable -type f | grep -q liftOver
then
    liftOver_exe=$additional_folder/liftOver
elif command -v liftOver &> /dev/null
then
    liftOver_exe=liftOver
else
    echo Please check that liftOver is installed.
    exit 1
fi

######## tabix

if find $additional_folder/ -executable -type f | grep -q tabix
then
    tabix_exe=$additional_folder/tabix
elif command -v tabix &> /dev/null
then
    tabix_exe=tabix
else
    echo Please check that tabix is installed.
    exit 1
fi

# DEFINE : $onekGP_directory

#### File creation 
echo 1000G_SNPs,1000G_SNPsDensity,aveMAF > data/1000g-freqsummary-$name.csv
echo gnomAD_SNP_count,gnomAD_SNP_density,gnomAD_avg_MAF > data/gnomad-$name.csv



###########################################################################################################################

# Obtain VCF Files from 1kGP 

###########################################################################################################################

#### Convert coordinates to hg19 genome version

awk -F',' 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1}' "$initial_data" > data/hg38-coordinates.bed     # Reformat coordinates for LiftOver

echo "$liftOver_exe data/hg38-coordinates.bed $chain_file data/hg19-coordinates.bed data/unlifted.bed &> /dev/null" >> errors.log
$liftOver_exe data/hg38-coordinates.bed $chain_file data/hg19-coordinates.bed data/unlifted.bed &> /dev/null

##

max="$last_rna_id"
var="$first_rna_id"

## If available, unzip previously VCF files, otherwise extract files
if [ -f data/1kGP-VCF.zip ]; then

    unzip data/VCF.zip &> /dev/null

else
    mkdir -p data/1kGP-VCF
    
    ## Obtain VCF files according to hg19 chromosome coordinates
    while [ $var -le $max ]
    do
        ## Grab the coordinates for each RNA as the counter increases
        line=$( grep -w "RNA$var" data/hg19-coordinates.bed | cut -f 1,2,3 | tr '\t' ' ' | tr -d "chr" | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' )
        
        ## If previous coordinate was associated with an annotated chromosome, use hg38 coordinates
        if [ -z $line ] || [[ $line == *"Un_"* ]]; then
      
            line=$( grep -w "RNA$var" $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | tr -d "chr" | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' )
            
        ## Reformatting scaffolds of annotated chromosomes
        elif [[ $line == *"_"* ]]; then
        
            a=$( echo $line | cut -d "_" -f 1 )
            b=$( echo $line | cut -d ":" -f 2 )
            line=$( echo $a:$b )

        fi
        
        rna_id='RNA'$var'-1k'
        chr=$( echo $line | cut -d ':' -f 1 )
	
        ## The X chromosome is downloaded from a different folder to the autosomal chromosomes.
        if [ $chr == 'X' ]; then
        
            #until
	        echo "$tabix_exe -f -h $onekGP_directory/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz $line > $rna_id.vcf 2>>tabix.log ;"  >> tabix.log		                                   
	        #$tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
	        $tabix_exe -f -h $onekGP_directory/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz $line > $rna_id.vcf 2>>tabix.log ;
	        #do
	        echo Processing $rna_id.vcf >> tabix.log ; echo >> tabix.log ;
	        #./data/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	        #./data/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	        #done
            ######## Tabix is very sensitive to crashing, so repeat the command until a valid VCF file has been downloaded.
        else
            #until
	        echo "$tabix_exe -f -h $onekGP_directory/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz $line > $rna_id.vcf 2>>tabix.log ;"  >> tabix.log
	        #$tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
	        $tabix_exe -f -h $onekGP_directory/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz $line > $rna_id.vcf 2>>tabix.log ;
	        #do
	        echo Processing $rna_id.vcf >> tabix.log ; echo >> tabix.log ;
	        #done
        fi                                                                            

        ######## Move VCF so it doesn't need to be redownloaded.
        mv $rna_id.vcf data/1kGP-VCF
        
        (( var++ ))
       
        rm -rf text.file
        
    done
fi

######## Remove excess files
rm -rf data/hg19-coordinated.bed
# rm -rf data/hg38-coordinates.bed
rm -rf data/unlifted.bed



############################################################################################################################

# Calculate features from 1kGP

############################################################################################################################

## Function declaration
VCF2summary() {

    local input_vcf="$1"
    local len="$2"

    ## Determine SNP frequencies for each ncRNA sequence
    echo "$vcftools_exe --vcf $input_vcf --freq --out data/vcftools-output &> /dev/null" >> tabix.log   
    $vcftools_exe --vcf $input_vcf --freq --out data/vcftools-output &> /dev/null

    ## Extract SNPs for each sequence at frequency observed at that position into a parsable format
    awk -F'\t' 'NR > 1 {print $5":"$6}' data/vcftools-output.frq > data/snps    

    SNPs_count=$( wc -l < data/snps)
    MAF_sum=$(awk -F':' '{sum += $4} END {print sum}' data/snps)

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
{
    read
    while IFS=, read -r id _ _ start end seq                              
    do
        len=$(( $end-$start ))                                              # to calculate SNP density
        input_vcf="data/1kGP-VCF/$id-1k.vcf"
        echo "[$id][$start][$end][$len][$id-1k.vcf]"  >> tabix.log  

        VCF2summary "$input_vcf" "$len"                                     

        echo "$SNPs_count,$SNPs_density,$MAF_average" >> data/1kGP-freqsummary-$name.csv
    
    done
   
} < $initial_data   

rm -rf data/hg19-coordinates.bed
rm -rf data/hg38-coordinates.bed
rm -rf data/snps
# remove VCF files or zip them? 

############################################################################################################################

# Obtain VCF Files from gnomAD and calculate features

############################################################################################################################

mkdir -p data/VCF-gnomad 

awk -F',' 'NR > 1 {print $3":"$4"-"$5}' "$initial_data" > data/coordinates

### Variables 
max="$last_rna_id"
var="$first_rna_id"

{
    read
    while IFS= read -r line
    do
        tabix_input="$line"

        start=$(echo "$line" | awk -F'[:-]' '{print $2}')
        end=$(echo "$line" | awk -F'[:-]' '{print $3}')                                 # For SNP density
        len=$(( $end-$start ))  

        rna_id="RNA$var-gnomad"

        chr=$(echo "$line" | awk -F'[:-]' '{print $1}')
        gnomad_database=data/raw/gnomad/gnomad.genomes.v4.0.sites.$chr.vcf.bgz          # make it a variable with user input

        echo Extracting $rna_id.vcf >> tabix.log ;
        echo "$tabix_exe -f -h $gnomad_database $tabix_input > $rna_id.vcf 2>>tabix.log ;" >>tabix.log
        $tabix_exe -f -h $gnomad_database $tabix_input > $rna_id.vcf 2>>tabix.log ;
        #do
        echo >> tabix.log 
#       done

        awk -F';' '{for(i=1;i<=NF;i++) if($i ~ /^AF=/) {split($i, a, "="); print a[2]}}' $rna_id.vcf > data/snps-gnomad
        
        SNPs_count=$( wc -l < data/snps-gnomad)
        MAF_sum=$(awk -F':' '{sum += $1} END {print sum}' data/snps-gnomad)

        if [[ "$SNPs_count" == "0" ]] || [ -z "$SNPs_count" ]; then               # If no SNPs were recorded, set everything to NA
    
            SNPs_count='NA'
            SNPs_density='NA'
            MAF_average='NA'
    
        else

            SNPs_density=$(echo "scale=4; $SNPs_count/$len" | bc )   
	        MAF_average=$(echo "scale=5; $MAF_sum/$SNPs_count" | bc )
    
        fi
   
        echo "$SNPs_count,$SNPs_density,$MAF_average" >> gnomad-$name.csv 
    
        (( var++ )) 

        mv "$rna_id.vcf" "data/VCF-gnomad/$rna_id.vcf"
    
    done
   
} < data/coordinates


######## Zip VCF files for further analyses
zip -r VCF VCF &> /dev/null
rm -rf VCF/