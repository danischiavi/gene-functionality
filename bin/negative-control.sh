
#### Variables and files
initial_data=$1

file_name=$(basename "${initial_data%.*}" | sed 's/coords-negative-control//')                  # Define name and directory for output file

left_negative_coords=data/"$file_name"left-negative-control.csv
right_negative_coords=data/"$file_name"right-negative-control.csv
negative_control=data/"$file_name"negative-control-dataset.csv

genome_csv=data/raw/GRCh38.p14_genome.csv
gencode_bed=data/raw/gencode-ncrna-annotation.bed 
uniprot_bed=data/raw/unipAliSwissprot.bed

bedtools_exe=bedtools


negative_control() {

    local distance_to_rna=$1

    chromo=$( echo $chr | tr -d "chr" )


    ## Generate null sequences upstream of sequence

    left_end=$(( "$start" - "$distance_to_rna" ))
    left_start=$((   "$left_end" - "$left_len" ))                                                               # length based on exons's length
    
    if [[ "$left_end" -lt 0 ]] || [[ "$left_start" -lt 0 ]]; then 
    
        left_sequence=                                                                                          # Blank to filter out if negative coord

    else 
    
        left_sequence=$( grep -w "chromosome $chromo" "$genome_csv" | cut -f 2 | cut -c$left_start-$left_end )  # Extract seq from genome file 
    
    fi
   
    if [ ! -z $left_sequence ] && [[ $left_sequence != *"N"* ]]; then                                           # If no sequence extracted, then remove.

        echo "$chr,$left_start,$left_end,$left_sequence" >> "$left_negative_coords"

    fi


    ## Generate null sequences downstream of sequence
    if [ ! -z "$right_len" ]; then                                            # Right_len will be empty for short-ncrna due to single exon

        right_start=$(( "$end" + "$distance_to_rna" ))
        right_end=$(( "$right_start" + "$right_len" ))
    
    else

        right_start=$(( "$end" + "$distance_to_rna" ))
        right_end=$(( "$right_start" + "$left_len" ))                        # left_len is the length of the only exon for short-ncrna

    fi

    if [[ "$right_start" -lt 0 ]] || [[ "$right_end" -lt 0 ]]; then 
    
        right_sequence= 
    
    else 
    
        right_sequence=$( grep -w "chromosome $chromo" $genome_csv | cut -f 2 | cut -c$right_start-$right_end ) 
    
    fi
    
    if [ ! -z $right_sequence ] && [[ $right_sequence != *"N"* ]]; then   

        echo "$chr,$right_start,$right_end,$right_sequence" >> "$right_negative_coords"

    fi

}


###########################################################################################################################

# Generate sequences upstream and downstream of functional gene

###########################################################################################################################

distances_to_rna=("4000000" "8000000" "12000000" "16000000" "20000000")

if [ ! -e "$right_negative_coords" ] || [ ! -e "$left_negative_coords" ]; then 

    while IFS=, read chr start end left_len right_len; do   #start exon 1, end last exon, length exon 2 and 3

        for position in "${distances_to_rna[@]}"; do negative_control "$position"; done

    done < "$initial_data"

fi


###########################################################################################################################

# Filter upstream (left) and downstream (right) negative control sequences using UniProt and GENCODE

###########################################################################################################################

filter_out_functional(){
    
    local side=$1

    ## Temporary files
    negative_coords=data/${file_name}${side}-negative-control.csv
    bed_coords=data/${file_name}${side}-coordinates.bed
    bed_intersect_uniprot=data/${file_name}${side}-intersect-uniprot.bed
    bed_intersect_gencode=data/${file_name}${side}-intersect-gencode.bed
    intersect_uniprot=data/${file_name}${side}-intersect-uniprot
    intersect_gencode=data/${file_name}${side}-intersect-gencode

    ######## Reformat data for bedtools
    grep -v "Start" "$negative_coords" | awk -F',' '{print $1 "\t" $2 "\t" $3}' > "$bed_coords"

    ######## Uniprot filtering
    "$bedtools_exe" intersect -a "$uniprot_bed" -b "$bed_coords" > "$bed_intersect_uniprot"

    ######## Extract chromosome coordinates of overlapping negative control sequences
    cut -f 1-3 "$bed_intersect_uniprot" | tr '\t' ',' > "$intersect_uniprot"
    #left_count=$( cat "$intersect" | wc -l )

    ######## GENCODE filtering
    "$bedtools_exe" intersect -a "$gencode_bed" -b "$bed_coords" > "$bed_intersect_gencode"

    ######## Extract chromosome coordinates of overlapping negative control sequences
    cut -f 1-3 "$bed_intersect_gencode" | tr '\t' ',' > "$intersect_gencode"

    ######## Remove negative control sequences that overlap with known functional sequences
    while IFS=',' read -r chr start end _; do

        coordinates="${chr},${start},${end}"

        if grep -qF "$coordinates" "$intersect_uniprot"; then
            :
        elif grep -qF "$coordinates" "$intersect_gencode"; then
            :
        else
            (( count++ ))
            echo "RNA$count,No,$chr,$start,$end" >> "$negative_control"
        fi

    done < "$negative_coords"

    rm -rf "$bed_coords"
    rm -rf "$bed_intersect_uniprot"
    rm -rf "$bed_intersect_gencode"
    rm -rf "$intersect_uniprot"
    rm -rf "$intersect_gencode"

}

if [ ! -e "$negative_control" ]; then

    echo "ID,Functional,Chromosome,Start,End,Sequence" > "$negative_control"
    count=$(wc -l < "$initial_data")                           # Start counting from last functional seq  
    # count=$(( $(wc -l < "$initial_data") - 1 )) # replace for this count when add header to the initial_data 
    filter_out_functional 'left'
    filter_out_functional 'right'

else

    echo "$negative_control already exist"

fi

######## Generate negative control FASTA file
#grep -v "Start" ./data/negative-control-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > ./data/negative-control-seq.fa

###########################################################################################################################

#rm -rf "$negative_coords"
