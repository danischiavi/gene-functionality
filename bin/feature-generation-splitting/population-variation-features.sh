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
additional_folder=$2 

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

## to alter count according to what number the RNA IDs start at
initial_max=$( tail -1 $initial_data | cut -d ',' -f 1 | tr -d RNA )
initial_var=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
initial_count=1

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


#### File creation 
echo 1000G_SNPs,1000G_SNPsDensity,aveMAF > 1000g-freqsummary.csv
echo gnomAD_SNP_count,gnomAD_SNP_density,gnomAD_avg_MAF > gnomad.csv


######## Convert coordinates to hg19 genome version
echo "$liftOver_exe coordinates-rnaid.bed $chain_file output.bed unlifted.bed &> /dev/null" >> errors.log
$liftOver_exe input.bed $chain_file output.bed unlifted.bed &> /dev/null


###########################################################################################################################

#### Obtain VCF Files from 1kGP 
max="$initial_max"
var="$initial_var"
count="$initial_count"

######## If available, unzip previously downloaded VCF files
if [ -f VCF.zip ]
then
    unzip VCF.zip &> /dev/null
else
    mkdir -p VCF
    
    ######## Downloaded VCF files according to hg19 chromosome coordinates
    while [ $var -le $max ]
    do
        ######## Grab the coordinates for each RNA as the counter increases
        line=$( grep -w "RNA$var" output.bed | cut -f 1,2,3 | tr '\t' ' ' | tr -d "chr" | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' )
        
        ######## If previous coordinate was associated with an annoated chromosome, use hg38 coordinates
        if [ -z $line ] || [[ $line == *"Un_"* ]]
        then
            line=$( grep -w "RNA$var" $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | tr -d "chr" | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' )
            
        ######## Reformatting scaffolds of annotated chromosomes
        elif [[ $line == *"_"* ]]
        then
            a=$( echo $line | cut -d "_" -f 1 )
            b=$( echo $line | cut -d ":" -f 2 )
            line=$( echo $a:$b )
        else
            :
        fi
        
        name='RNA'$var'-1k'
        chr=$( echo $line | cut -d ':' -f 1 )
	
        ######## The X chromosome is downloaded from a different folder to the autosomal chromosomes.
        if [ $chr == 'X' ]
        then
            #until
	    echo "$tabix_exe -f -h $additional_folder/1000genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ;"  >> tabix.log		                                   
	    #$tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
	    $tabix_exe -f -h $additional_folder/1000genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ;
	    #do
	    echo Processing $name.vcf >> tabix.log ; echo >> tabix.log ;
	    #./data/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	    #./data/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	    #done
            ######## Tabix is very sensitive to crashing, so repeat the command until a valid VCF file has been downloaded.
        else
            #until
	    echo "$tabix_exe -f -h $additional_folder/1000genomes/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ;"  >> tabix.log
	    #$tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
	    $tabix_exe -f -h $additional_folder/1000genomes/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ;
	    #do
	    echo Processing $name.vcf >> tabix.log ; echo >> tabix.log ;
	    #done
        fi                                                                            

        ######## Move VCF so it doesn't need to be redownloaded.
        mv $name.vcf VCF/
        var=$((var+1))
        rm -rf text.file
        count=$(( $count + 1 ))
    done
fi

######## Remove excess files
rm -rf output.bed
rm -rf input.bed
rm -rf unlifted.bed



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

############################################################################################################################

# Calculate features from 1kGP

############################################################################################################################

{
    read
    while IFS=, read -r id _ _ start end seq        # Loop through VCF files
    do
        len=$(( $end-$start ))   # Sequence length
        f="$id-1k.vcf"
        input_vcf="VCF/$f"
        output_prefix="freq-afr"
        echo "[$id][$start][$end][$len][$f]"  >> tabix.log  

        VCF2summary "$vcftools_exe" "$input_vcf" "$output_prefix" "$len"

        echo "$count,$density,$average" >> 1000g-freqsummary.csv
    
    done
   
} < $initial_data   #would be better to loop through coordinates file? 



############################################################################################################################

# Obtain VCF Files from gnomAD and calculate features

############################################################################################################################

### Variables 
max="$initial_max"
var="$initial_var"
count="$intial_count"

{
    read
    while IFS='/t' read -r id start end
    do
        len=$(( $end-$start ))   # Sequence length
        name='RNA'$var'-gnomad'
        tabix_input=$id':'$start'-'$end

        echo Extracting $name.vcf >> tabix.log ;
        echo "$tabix_exe -f -h $gnomad_database $tabix_input > $name.vcf 2>>tabix.log ;" >>tabix.log
        $tabix_exe -f -h $gnomad_database $tabix_input > $name.vcf 2>>tabix.log ;
        #do
        echo >> tabix.log ;
#       done
    
        count=$( grep -v "#" $name.vcf | wc -l )
        density=$(echo "scale=3; $count/$len" | bc)
        maf=0
        
    for snp in $( grep -v "#" $name.vcf )  #  careful with IFS! maybe is worth to define it for this specific loop
    do
        af=$( echo $snp | cut -f 8 | cut -d ';' -f 3 | tr -d "AF=" ) >> af.txt #bcftools??
        af_decimal=$(printf "%.10f" "$af")                              #from scientific annotation to decimal for bc (check smallest maf in database)
        maf=$(echo "$maf+$af_decimal" | bc )
    done

    [[ "$count" -eq 0 ]] && avg_maf=NA || avg_maf=$(echo "scale=3; $maf/$count" | bc)

    [[ "$count" -eq 0 ]] && echo NA,NA,NA >> $date-gnomad.csv || echo $count,$density,$avg_maf >> gnomad.csv

    #rm -rf $name.vcf
    var=$((var+1))

    mv $name.vcf VCF/
    
    done
   
} < $coordinates