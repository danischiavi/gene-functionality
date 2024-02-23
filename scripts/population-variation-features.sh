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
chain_file=$3
onekGP_directory=$4
gnomad_directory=$5

first_rna_id=$(awk -F',' 'NR==2 {print $1}' "$initial_data" | tr -d 'RNA')
last_rna_id=$(awk -F',' 'END {print $1}' "$initial_data" | tr -d 'RNA') 

output_directory=data/population
mkdir -p "$output_directory"

# Final output file
file_name=${output_directory}/$(basename "${initial_data%.*}" | sed 's/-dataset//')
output_file_variation="$file_name"-variation.csv 

# Temporary output file
output_file_1kGP="$file_name"-1kGP-variation.csv                  
output_file_gnomAD="$file_name"-gnomAD-variation.csv 
 
# Variables
liftOver_exe=liftOver
tabix_exe=tabix
vcftools_exe=vcftools 


###########################################################################################################################

# Function declaration

###########################################################################################################################

VCF2summary() {

    local input_vcf="$1"
    local len="$2"

    ## Determine SNP frequencies for each ncRNA sequence
    output_vcftools="$file_name"-vcftools-output
    echo "$vcftools_exe --vcf $input_vcf --freq --out $output_vcftools &> /dev/null" >> tabix.log   
    $vcftools_exe --vcf "$input_vcf" --freq --out "$output_vcftools" &> /dev/null

    ## Extract SNPs for each sequence at frequency observed at that position into a parsable format
    snps_file="$file_name"-snps
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

###########################################################################################################################

# Extract features 

###########################################################################################################################
VCF_directory_1kGP="$file_name"-1kGP-VCF 

if [ ! -s "$output_file_1kGP" ]; then

    #### Obtain VCF Files from 1kGP ####

    # If available, unzip previously VCF files, otherwise extract files   
    if [ ! -e "$VCF_directory_1kGP" ]; then
    
        if [ -f "$VCF_directory_1kGP".zip ]; then

            unzip "$VCF_directory_1kGP".zip &> /dev/null

        else

            mkdir -p "$VCF_directory_1kGP"

            # Convert coordinates to hg19 genome version
            current_coord="$file_name"-hg38-coordinates.bed
            old_coord="$file_name"-hg19-coordinates.bed
            unlifted_coord="$file_name"-unlifted.bed

            # Reformat coordinates for LiftOver
            awk -F',' 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1}' "$initial_data" > "$current_coord"  

            echo "$liftOver_exe $current_coord $chain_file $old_coord $unlifted_coord &> /dev/null" >> errors.log
            "$liftOver_exe" "$current_coord" "$chain_file" "$old_coord" "$unlifted_coord" &> /dev/null

            max="$last_rna_id"
            var="$first_rna_id"
    
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
        
                rna_id="$VCF_directory_1kGP"/'RNA'"$var"
                chr=$( echo $line | cut -d ':' -f 1 )
	
                ## The X chromosome is downloaded from a different folder to the autosomal chromosomes.
                if [ $chr == 'X' ]; then
        
                    #until
	                echo "$tabix_exe -f -h $onekGP_directory/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz $line > $rna_id.vcf 2>>tabix.log ;"  >> tabix.log		                                   
	                #$tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
	                "$tabix_exe" -f -h "$onekGP_directory"/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz "$line" > "$rna_id".vcf 2>>tabix.log ;
	                #do
	                echo Processing $rna_id.vcf >> tabix.log ; echo >> tabix.log ;
	                #./data/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	                #./data/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	                #done
                    ######## Tabix is very sensitive to crashing, so repeat the command until a valid VCF file has been downloaded.

                else
                    #until
	                echo "$tabix_exe -f -h $onekGP_directory/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz "$line" > "$rna_id".vcf 2>>tabix.log ;"  >> tabix.log
	                #$tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
	                "$tabix_exe" -f -h "$onekGP_directory"/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz "$line" > "$rna_id".vcf 2>>tabix.log ;
	                #do
	                echo Processing "$rna_id".vcf >> tabix.log ; echo >> tabix.log ;
	                #done
                fi                                                                            
        
                (( var++ ))
       
                rm -rf text.file
        
            done

        fi
    fi
    #### Calculate features from 1kGP VCF files ####

    echo "1kGP_SNPs,1kGP_SNP_density,1kGP_aveMAF" > "$output_file_1kGP"

    tail -n +2 "$initial_data"  | while IFS=, read -r id _ _ start end seq; do                              
    
        len=$(( $end-$start ))                                                  # to calculate SNP density
        input_vcf="$VCF_directory_1kGP"/"$id".vcf
        echo "[$id][$start][$end][$len]["$id".vcf]"  >> tabix.log  

        VCF2summary "$input_vcf" "$len"                                         # Call function declared above           

        echo "$SNPs_count,$SNPs_density,$MAF_average" >> "$output_file_1kGP"
    
    done

    ## Zip VCF files for further analyses
    zip -r "$VCF_directory_1kGP".zip "$VCF_directory_1kGP" &> /dev/null  && rm -rf "$VCF_directory_1kGP"

fi


############################################################################################################################

# Obtain VCF Files from gnomAD and calculate features

############################################################################################################################
VCF_directory_gnomAD="$file_name"-gnomAD-VCF

if [ ! -s "$output_file_gnomAD" ]; then
    
    echo "gnomAD_SNPs,gnomAD_SNP_density,gnomAD_aveMAF" > "$output_file_gnomAD"
    
    if [ ! -e "$VCF_directory_gnomAD" ]; then
    
        if [ -f "$VCF_directory_gnomAD".zip ]; then

            unzip "$VCF_directory_gnomAD".zip &> /dev/null
        
        else

            mkdir -p "$VCF_directory_gnomAD"

        fi

    fi

    coordinates="$file_name"-tabix-coordinates

    awk -F',' 'NR > 1 {print $3":"$4"-"$5}' "$initial_data" > "$coordinates"

    ### Variables 
    max="$last_rna_id"
    var="$first_rna_id"

    while IFS= read -r line; do
        
        rna_id="$VCF_directory_gnomAD"/RNA"$var"
        
        if [ ! -e "$rna_id" ]; then

            tabix_input="$line"

            start=$(echo "$line" | awk -F'[:-]' '{print $2}')
            end=$(echo "$line" | awk -F'[:-]' '{print $3}')                                         
            len=$(( $end-$start ))                                                                     # For SNP density 

            chr=$(echo "$line" | awk -F'[:-]' '{print $1}')
            gnomad_database="$gnomad_directory"/gnomad.genomes.v4.0.sites."$chr".vcf.bgz               # make it a variable with user input

            echo Extracting :"$rna_id".vcf >> tabix.log ;
            echo "$tabix_exe -f -h $gnomad_database $tabix_input > $rna_id.vcf 2>>tabix.log ;" >>tabix.log
            $tabix_exe -f -h "$gnomad_database" "$tabix_input" > "$rna_id".vcf 2>>tabix.log ;
            #do
            echo >> tabix.log 
#           done
        fi

        snps_file_gnomAD="$file_name"-snps-gnomAD
        awk -F';' '{for(i=1;i<=NF;i++) if($i ~ /^AF=/) {split($i, a, "="); print a[2]}}' "$rna_id".vcf > "$snps_file_gnomAD"
        
        SNPs_count=$( wc -l < "$snps_file_gnomAD")
        MAF_sum=$(awk -F':' '{sum += $1} END {printf "%.10f", sum}' "$snps_file_gnomAD")
      

        if [[ "$SNPs_count" == "0" ]] || [ -z "$SNPs_count" ]; then                                     # If no SNPs were recorded, set everything to NA
    
            SNPs_count='NA'
            SNPs_density='NA'
            MAF_average='NA'
    
        else

            SNPs_density=$(echo "scale=6; $SNPs_count/$len" | bc )   
	        MAF_average=$(echo "scale=6; $MAF_sum/$SNPs_count" | bc )
    
        fi

        
        echo "$SNPs_count,$SNPs_density,$MAF_average" >> "$output_file_gnomAD" 
    
        (( var++ )) 
    
    done < "$coordinates"

    ## Zip VCF files for further analyses
    zip -r "$VCF_directory_gnomAD".zip "$VCF_directory_gnomAD" &> /dev/null && rm -rf "$VCF_directory_gnomAD" 

fi

if [ ! -s "$output_file_variation" ]; then

    paste -d',' "$output_file_1kGP" "$output_file_gnomAD"   > "$output_file_variation"

fi
 

### Remove excess files
rm -rf "$output_vcftools".frq
rm -rf "$current_coord"
rm -rf "$old_coord"
rm -rf "$unlifted_coord"
rm -rf "$snps_file" 
rm -rf "$output_vcftools"
rm -rf "$coordinates"
rm -rf "$snps_file_gnomAD"
rm -rf "$snps_file" 