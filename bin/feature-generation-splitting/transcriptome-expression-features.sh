#!/bin/bash 

# Script Name: transcriptome-expression-features.sh
#
# Description: This script calculates Maximum read depth (MRD) and the average number of reads per kilobase of transcripts per million mapped reads (RPKM) for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences.
#
# Input: $1 is the file identifier for the dataset and fasta file. Eg: 200702-functional-ncrna
#        $2 is the location of the folder with the local databases/datasets or version specific executables
# 
# Any additional required files, directories or dependencies will be requested when the script is run and require manual
# input.
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.


############################################################################################################################

initial_data=$1

output_directory=data/transcriptome
mkdir -p "$output_directory"

# Variables
encode_folder=/Volumes/archive/userdata/student_users/danielaschiavinato/dani-scratch/features-of-function-data
samtools_exe=samtools

#### Function declaration 
run_encode() {

    local sample_type=$1

    output_file="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')"$sample_type".csv"   # Define name and directory for output file
    echo RPKM_$sample_type,MRD_$sample_type > $output_file

    ######## Define location of downloaded ENCODE RNAseq data
    encode_data=$encode_folder/ENCODE-$sample_type/*.bam
    total_read_file=$encode_folder/ENCODE-$sample_type/total-read-count.txt 

    {
        read
        while IFS=',' read -r _ _ chr start end _; do
    
            length=$(( $end-$start ))
            length_kb=$( echo "scale=2; $length/1000" | bc )                                             # Convert to Kb
            input_samtools=$chr':'$start'-'$end
    
        ######## Record RPKM and MRD across all RNAseq files
            count_max=0                                                                                  # To record maximum RPKM and MRD observed
            max_depth=0                                                                                  # To record maximum MRD observed

            for file in $encode_data; do                                                                     # Loop through all RNAseq datasets
     
                read=$( $samtools_exe view -c $file $input_samtools 2>>errors.log )                      # Number of reads for sequence of interest
                file_name=$( echo $file | rev | cut -d '/' -f 1 | rev )                                  # ID for RNAseq dataset
                total_read=$( grep "$file_name" $total_read_file | cut -d ',' -f 2 )                     # Obtain total number of reads for specific RNAseq dataset (precomputed to save time)
                scaling=$(   echo "$total_read/1000000" | bc )                                           # Create the "by per million" scaling 
                rpm=$(   echo "scale=2; $read/$scaling" | bc )                                           # Reads per million
                rpkm=$( echo "scale=2; $rpm/$length_kb" | bc )                                           # Final RPKM calculation
        
                # If RPKM is greater than previously recorded, then overwrite the variable.
                if (( $(echo "$rpkm > $count_max" | bc -l) )) ; then count_max=$rpkm ; else : ; fi
        
                max=$( $samtools_exe depth -r $input_samtools $file 2>>errors.log | cut -f 3 | sort -V | tail -1 )  # Maximum Read Depth
       
                # If MRD is greater than previously recorded, then overwrite the variable.
                if [[ $max -ge $max_depth ]] ; then max_depth=$max ; fi

            done 

            [ -z "$max_depth" ] && max_depth=0      
                                                                  # Record MRD as zero rather than blank
            if [ -z "$count_max" ] ; then count_max=NA ; max_depth=NA ; fi                       # But if no transcription whatsoever, set everything to NA 

            echo $count_max,$max_depth >> $output_file

        done 

    } < $initial_data

}

############################################################################################################################

# Extracting expression for a region across Total & Small ENCODE datasets

############################################################################################################################

if [ ! -e ${output_directory}/*-tissue.csv ]; then run_encode 'tissue'; fi

if [ ! -e ${output_directory}/*-primary-cell.csv ]; then run_encode 'primary-cell'; fi

paste -d',' ${output_directory}/*-tissue.csv ${output_directory}/*-primary-cell.csv > "${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')transcriptome.csv"  