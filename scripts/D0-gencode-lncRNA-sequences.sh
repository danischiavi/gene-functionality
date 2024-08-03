#!/bin/bash
#
# Script Name: non-coding-sequences.sh
#
# Author: Daniela Schiavinato
# Last edited: April 2024
#
# Description: 

#   - Generate Functional datasets: 
#                                 functional-lncrna-exon1 & functional-lncrna-exon2
#
#   - Extracts coordinates to generate negative control dataset downstream using negative-control.sh:                     
#                                 lncrna-coords-negative-control

###########################################################################################################################
#### GENERAL SETUP ####
###########################################################################################################################

##### Input files ##### 
gencode_annotations=$1
genome_seq=$2

# Gene list 
awk -F'\t' '$3 == "gene" {split($NF, arr, ";"); for (i in arr) {if (match(arr[i], /^ID=([^;]+)/, id)) {print id[1]}}}' "$gencode_annotations" > data/gencode-lncRNA/gencode-genes-list 	# maybe add it to readmedata 

#### Final Output files #####
dir='data/datasets'
mkdir -p "$dir"

# Functional datasets 
exon_one="${dir}/gencode-lncrna-exon1-dataset.csv"
exon_two="${dir}/gencode-lncrna-exon2-dataset.csv" 

#### Temporary files ####
exon_one_info="${dir}/gencode-lncrna-exon1-info"
exon_two_info="${dir}/gencode-lncrna-exon2-info"
exon_one_tmp="${dir}/gencode-lncrna-exon1-tmp"
exon_two_tmp="${dir}/gencode-lncrna-exon2-tmp"

exon_one_unsorted="${dir}/gencode-lncrna-exon1-unsorted.csv" 
exon_two_unsorted="${dir}/gencode-lncrna-exon2-unsorted.csv"  

##### Constrains ##### 
lower_limit='74'              # Given by the size distribution analysis (10&90% percentile)
upper_limit='1500'

###########################################################################################################################
#### FUNCTION DECLARATION ####
###########################################################################################################################
# Get info from gff3 file 
gff2Info() {                                                                                

    local exons=$1
	local id=$2

    coords_one=$( awk 'NR==1 {print $1, $4, $5, $7}' "$exons" )                                                  
    coords_two=$( awk 'NR==2 {print $1, $4, $5}' "$exons" )                                                    
      												
    IFS=' ' read -r chr start_one end_one strand <<< "$coords_one"                                               
    IFS=' ' read -r _ start_two end_two <<< "$coords_two"
	start_one_bed=$(( start_one - 1 ))															
	start_two_bed=$(( start_two - 1 ))
                                                                       
    len_one=$(( end_one - start_one_bed ))
	len_two=$(( end_two - start_two_bed ))           
       
	# Exclude sequences out of length limits 
	if { [ "$len_two" -ge "$lower_limit" ] && [ "$len_two" -le "$upper_limit" ]; } && \
    { [ "$len_one" -ge "$lower_limit" ] && [ "$len_one" -le "$upper_limit" ]; }; then 
                                        
		echo -e "${chr}\t${start_one}\t${end_one}\t${id}\t.\t${strand}" >> "$exon_one_info"		# following bed format 
        echo -e "${chr}\t${start_two}\t${end_two}\t${id}\t.\t${strand}" >> "$exon_two_info"        
    fi 
}


reformat_file() {

    local input_file=$1
    local output_file=$2

    awk -F'\t' '
    function calc_ambiguity(seq,   i, amb_nucleotides, total_nucleotides, ambiguity_percent) {
        amb_nucleotides = 0
        total_nucleotides = length(seq)
        
        for (i = 1; i <= total_nucleotides; i++) {
            if (index("RYMKSWBDHVNrymkswbdhvn", substr(seq, i, 1))) {
                amb_nucleotides++
            }
        }
        
        ambiguity_percent = (amb_nucleotides * 100) / total_nucleotides
        return ambiguity_percent
    	
		}

    {
        split($1, parts, "::")
        split(parts[2], coords, ":")
        gene_id = parts[1]
        chromosome = coords[1]
        split(coords[2], range, "-")
        start = range[1] + 1 
        end = substr(range[2], 1, index(range[2], "(") -1)

        sequence = $2
        ambiguity = calc_ambiguity(sequence)

        if (ambiguity <= 5) {
            printf("NA,%s,%s,%s,%s,%s\n", chromosome, start, end, sequence, gene_id)
        }
    }' "$input_file" >> "$output_file"

	# start + 1 since is following bed format 

    rm -rf "$input_file"
}

###########################

sort_output(){

    local input_file=$1
    local output_file=$2

	local file_id=data/datasets/"$(basename "${input_file%.*}" | sed 's/-unsorted.csv//')" 

	sort -t ',' -k2,2 -k3,3n -k4,4n "$input_file" > "$file_id"-sorted-columns

	awk -v start=1 -v end=$(( $(wc -l < "$input_file")))  '
    BEGIN {
        for (i = start; i <= end; i++) {
            print "RNA" i
        }
    }' > "$file_id"-id-column

    (echo "ID,Functional,Chromosome,Start,End,Sequence,GeneID"; paste -d ',' "$file_id"-id-column "$file_id"-sorted-columns) > "$output_file"

	rm -rf "$file_id"-id-column "$file_id"-sorted-columns
	rm -rf "$input_file"
}


###########################################################################################################################
#### EXTRACT SEQUENCES ####
###########################################################################################################################

if [ ! -s "$exon_one" ] || [ ! -s "$exon_two" ]; then

	if [ ! -s "$exon_one_info" ] || [ ! -s "$exon_two_info" ]; then  
    
		declare -a "IDs_lncrna=()"
    	mapfile -t IDs_lncrna < data/gencode-genes-list

		for lncrna_id in "${IDs_lncrna[@]}"; do

			awk -F'\t' -v id="$lncrna_id" '$0 ~ "gene_id="id && $3 == "exon"' "$gencode_annotations" > data/gencode-lncRNA/gencode-exons

            if [ "$(wc -l < data/gencode-exons)" -ge 2 ]; then gff2Info data/gencode-lncRNA/gencode-exons "$lncrna_id"; fi    
		
    	done
	fi

	
	## Extract sequences from genome ##
	if [ ! -s "$exon_one_tmp" ]; then  
		bedtools getfasta -fi "$genome_seq" -bed "$exon_one_info" -fo "$exon_one_tmp" -s -name -tab  
			# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
			# -name Use the name field and coordinates for the FASTA header
			# tab Report extract sequences in a tab-delimited format instead of in FASTA format.
	fi 

	if [ ! -s "$lncrna_exon_two_tmp" ]; then 
		bedtools getfasta -fi "$genome_seq" -bed "$exon_two_info" -fo "$exon_two_tmp" -s -name -tab 
	fi

	## Format Final file ## 
	if [ ! -s "$exon_one_unsorted" ]; then 
		reformat_file "$exon_one_tmp" "$exon_one_unsorted"
		sort_output "$exon_one_unsorted" "$exon_one"
	fi

	if [ ! -s "$exon_two_unsorted" ]; then 
		reformat_file "$exon_two_tmp" "$exon_two_unsorted"
		sort_output "$exon_two_unsorted" "$exon_two"
	fi

else 
	echo "Files already exist" 
fi 
    

rm -rf data/gencode-genes-list 
rm -rf data/gencode-exons