#!/bin/bash
# 
## This script generates the corresponding negative control dataset for the input file (protein-codingRNA, lncRNA or short-ncRNA)       
# - Outputs one file for each exon (Exon 2&3 for protein-coding; 1&2 for lncRNA and single exon for short-ncRNA)      
# - Note that A refers to the exon of a functional seq of the RNA type (Exon 1 for lncRNA and Exon 2 for protein-codingRNA);
# B to the second exon (Exon 2 for lncRNA and Exon 3 for protein-codingRNA); and
# single refers to sequences which have a single exon (short-ncRNA)
#  
######## General Idea ########
# 1) Reading line by line, define a sampling region from the gene to a fix value (first downstream) 
#    This fix value corresponds to the furthest away we want to sample a negative control (highest d - see line 15). 
# 2) Find the intersection between the sampling region and the "genes complement"
# 3) For each exon, find the closest region of the "genes complement" to a position=(Exon + d )  
# 4) Take a random sample within this closest region of length=(exon's length)
# 6) Repeat (3) & (4) for each d. Note values for d are pre-defined - see manuscript for justification 
# 7) Repeat (1)-(6) for upstream 
##############################

#### General Set Up ####

## Input files ##
# Of each genes from were functional seq was selected: chr, start coords of first exon, last coord of last exon, length exonA and length exonB  
initial_data=$1  
       
# csv file with genome sequences for each chromosome 
genome_seq=$2
# regions of genome lacking of genes (annotated from RNAcentral and Gencode)
genes_complement=$3

initial_data_seq_numb=$(( $(wc -l < "$initial_data") ))    
# distances from genome to find a region to sample the negative control. See manuscript for justification
distances_to_seq=("1000" "10000" "100000" "1000000" "5000000")                             

# Name for files
file_name=data/datasets/$(basename "${initial_data%.*}" | sed 's/-coords-negative-control//')                 

## Temporary files ##
A_info="$file_name"-A-info
B_info="$file_name"-B-info
A_tmp="$file_name"-A-tmp
B_tmp="$file_name"-B-tmp
A_unsorted="$file_name"-A-unsorted
B_unsorted="$file_name"-B-unsorted
short_info="$file_name"-info
short_tmp="$file_name"-tmp
short_unsorted="$file_name"-unsorted

sampling_region="$file_name"-sampling-region
complement_sampling_region="$file_name"-complement-sampling-region
                                                                                                                                   
closest_input="$file_name"-closest-input
closest_output="$file_name"-closest-output
random_output="$file_name"-random-coords

## Final output files ##
# Note the output files depends on the initial_data input file 
negative_control_single="$file_name"-negative-control-dataset.csv   # File for short-ncrna is the final negative control dataset
negative_control_A="$file_name"-exonA-negative-control.csv          # File for lncrna and protein will be rename after this script to associated with the corresponding exon  
negative_control_B="$file_name"-exonB-negative-control.csv

###########################################################################################################################
#### FUNCTION DECLARATION ####
###########################################################################################################################
# Finds the intersection between the sampling region and the "genes complement", removing the regions which length is smaller than desired 
# $genes_complement is a global variable of the script which is assigned with a script input

intersect() {

	sampling_region=$1
	len=$2
	complement_sampling_region=$3

	bedtools intersect -a "$sampling_region" -b "$genes_complement" \
	| awk -F'\t' -v len="$len" '{ if (($3-$2) > len) print $0"\t"$3-$2}' \
	| sort -k1,1 -k2,2n -k3,3n \
	> "$complement_sampling_region"

}


###########################
# Takes a random sample within a region:  
# selects a random start and end coordinate for the negative control from a specific region (closest to position=(exon+d))
sample_random() {
    
	local chrRandom="$1"
    local startRandom="$2"
    local endRandom="$3"
    local lenRandom="$4"

	rand_end=$(( endRandom + 1 ))

	while [ "$rand_end" -gt "$endRandom" ]; do 
    
		# Generate random start coordinate within the specified range
    	rand_start=$(( RANDOM % (endRandom - startRandom - lenRandom + 2) + startRandom ))      # 2 accounts for start being inclusive 

    	# Calculate end coordinate based on start and region length
    	rand_end=$(( rand_start + lenRandom - 1 ))

	done

    # Print sampled region
    echo -e "$chrRandom\t$rand_start\t$rand_end" > "$random_output" 
}

###########################
# Organize auxiliary functions and defines the positions=(exon+d)
#### For sequences with more than 1 exon #### 
multiple_exons_negative_control() {

    flow=$1
    sampling_regions=$2

    len_region=$(( lenA + lenB )) 

    for d in "${distances_to_seq[@]}"; do    
    
        # Downstream
        if [ "$flow" == 'downstream' ]; then
            
            initial_end=$(( end_gene + d ))
            initial_start=$(( initial_end - len_region ))
        
        # Upstream
        elif [ "$flow" == 'upstream' ]; then 

            initial_start=$(( start_gene - d ))
            if [ "$initial_start" -lt 0 ]; then initial_start=0; fi             # To stay within the chromosome region (downstream is limitated by the intersect of the complement)

            initial_end=$(( initial_start + len_region ))
        
        fi

    if [ "$initial_end" -gt 0 ] && [ "$initial_start" -gt 0 ]; then                         # Within chromosome limits  
            echo -e "$chr\t$initial_start\t$initial_end" > "$closest_input" 

            # Find closest region of genes complement 
            bedtools closest -a "$closest_input" -b "$sampling_regions" -D ref -t first  > "$closest_output"
        
            # Find random coordinates within the closest region
            IFS=$'\t' read -r _ _ _ chrClosest startClosest endClosest _ < "$closest_output" 

            #### EXON A ####
            sample_random "$chrClosest" "$startClosest" "$endClosest" "$lenA" 
            
            # Redefine variables with random coordinates 
            IFS=$'\t' read -r chrClosestRdm startClosestRdm endClosestRdm < "$random_output" 

			# Calculate distance to gene border
   			if [ "$flow" == 'downstream' ]; then distance_gene_A=$(( startClosestRdm - end_gene ));
   			elif [ "$flow" == 'upstream' ]; then distance_gene_A=$(( start_gene - endClosestRdm ));
   			fi

			(( countA++ ))
			echo -e "$chrClosestRdm\t$startClosestRdm\t$endClosestRdm\tRNA$countA.$distance_gene_A\t.\t$strand" >> "$A_info"

            #### EXON B ####
            sample_random "$chrClosest" "$startClosest" "$endClosest" "$lenB" 
            
            # Redefine variables with random coordinates 
            IFS=$'\t' read -r chrClosestRdm startClosestRdm endClosestRdm < "$random_output" 

			# Calculate distance to gene border
   			if [ "$flow" == 'downstream' ]; then distance_gene_B=$(( startClosestRdm - end_gene ));
   			elif [ "$flow" == 'upstream' ]; then distance_gene_B=$(( start_gene - endClosestRdm ));
   			fi

			(( countB++ ))
            echo -e "$chrClosestRdm\t$startClosestRdm\t$endClosestRdm\tRNA$countB.$distance_gene_B\t.\t$strand" >> "$B_info"

        fi
        
    done

    rm -rf "$closest_input"
    rm -rf "$closest_output"
}
    
#### For sequences with single exon (short-ncRNA) #### 
single_exon_negative_control() {

    flow=$1
    sampling_regions=$2
   
    for d in "${distances_to_seq[@]}"; do                                                       # Defined on line 15
    
        # Downstream
        if [ "$flow" == 'downstream' ]; then
            
            initial_end=$(( end_gene + d ))
            initial_start=$(( initial_end - len ))
        
        # Upstream
        elif [ "$flow" == 'upstream' ]; then 

            initial_start=$(( start_gene - d ))
            if [ "$initial_start" -lt 0 ]; then initial_start=0; fi 

            initial_end=$(( initial_start + len ))
        
        fi

        if [ "$initial_end" -gt 0 ] && [ "$initial_start" -gt 0 ]; then                         # Within chromosome limits  
            echo -e "$chr\t$initial_start\t$initial_end" > "$closest_input" 

            # Find closest region of genes complement 
            bedtools closest -a "$closest_input" -b "$sampling_regions" -D ref -t first  > "$closest_output"
			
			IFS=$'\t' read -r _ _ _ chrClosest startClosest endClosest _ < "$closest_output"
            
			### xxxxxxx
             
            sample_random "$chrClosest" "$startClosest" "$endClosest" "$len" 

            # Redefine variables with random coordinates 
            IFS=$'\t' read -r chrClosestRdm startClosestRdm endClosestRdm < "$random_output" 

			if [ "$flow" == 'downstream' ]; then distance_gene_single=$(( startClosestRdm - end_gene ));
   			elif [ "$flow" == 'upstream' ]; then distance_gene_single=$(( start_gene - endClosestRdm ));
   			fi

			(( countSingle++ ))
			echo -e "$chrClosestRdm\t$startClosestRdm\t$endClosestRdm\tRNA$countSingle.$distance_gene_single\t.\t$strand" >> "$short_info"
        fi
        
    done

    rm -rf "$closest_input"
    rm -rf "$closest_output"
}
    
###########################
## Format Final file ## 
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
		split(parts[1], header_parts, ".")
        distance = header_parts[2]
        split(parts[2], coords, ":")
        chromosome = coords[1]
        split(coords[2], range, "-")
        start = range[1] + 1 
        end = substr(range[2], 1, index(range[2], "(") - 1)

        sequence = $2
        ambiguity = calc_ambiguity(sequence)

        if (ambiguity <= 5) {
            printf("No,%s,%s,%s,%s,%s\n", chromosome, start, end, sequence, distance)
        }
    }' "$input_file" >> "$output_file"

    rm -rf "$input_file"
}


###########################
# Remove duplicates, sort sequences and rename (RNAid) 

sort_output(){

	local input_file=$1
    local output_file=$2

    local file_id=data/datasets/"$(basename "${input_file%.*}" | sed 's/-unsorted//')" 

	sort -t ',' -k2,2 -k3,3n -k4,4n "$input_file" > "$file_id"-sorted-columns

    numb_seqs=$(wc -l < "$file_id"-sorted-columns)

	awk -v start="$initial_data_seq_numb" -v end="$(( initial_data_seq_numb + numb_seqs - 1 ))" '
    BEGIN {
        for (i = start; i <= end; i++) {
            print "RNA" i
        }
    }' > "$file_id"-id-column

    (echo "ID,Functional,Chromosome,Start,End,Sequence,DistanceGene"; paste -d ',' "$file_id"-id-column "$file_id"-sorted-columns) > "$output_file"

	rm -rf "$file_id"-id-column "$file_id"-sorted-columns 
	
}

###########################################################################################################################
#### EXTRACT NEGATIVE CONTROLS ####
###########################################################################################################################
 
# Define if it's a single or multiple exon dataset counting the fields of initial data file (4: short-ncRNA; 5: lncRNA and protein-codingRNA)
num_fields=$( awk -F',' '{print NF}' "$initial_data" | sort -nu | tail -n 1 )                      

#### Single exon sequences ####
if [[ "$num_fields" -eq 5 ]]; then                                                                      

    if [ ! -s "$negative_control_single" ]; then 

        countSingle=$(( $(wc -l < "$initial_data") - 1 ))                                           # Starts counting from last corresponding functional seq  
        
        tail -n +2 "$initial_data"  | while IFS=, read -r chr start_gene end_gene len strand; do   
          
            # Downstream 
            end_sampling_region=$(( end_gene + 5000000 ))
            echo -e "$chr\t$end_gene\t$end_sampling_region" > "$sampling_region"

            intersect "$sampling_region" "$len" "$complement_sampling_region"
           
            single_exon_negative_control 'downstream' "$complement_sampling_region"
            
            # Upstream 
            start_sampling_region=$(( start_gene - 5000000 ))
            if [ "$start_sampling_region" -lt 0 ]; then start_sampling_region=0; fi 

            echo -e "$chr\t$start_sampling_region\t$start_gene" > "$sampling_region"

            intersect "$sampling_region" "$len" "$complement_sampling_region"
           
            single_exon_negative_control 'upstream' "$complement_sampling_region" 

			rm -rf "$complement_sampling_region" "$sampling_region" "$random_output"

        done 

		bedtools getfasta -fi "$genome_seq" -bed "$short_info" -fo "$short_tmp" -s -name -tab 
		# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
		# -name Use the name field and coordinates for the FASTA header
		# tab Report extract sequences in a tab-delimited format instead of in FASTA format

		reformat_file "$short_tmp" "$short_unsorted"
		
		sort_output "$short_unsorted" "$negative_control_single"  

		rm -rf "$short_info" "$short_tmp" "$short_unsorted" 
    
	fi

fi


# Multiple exons sequences #
if [[ "$num_fields" -eq 6 ]]; then 

    if [ ! -s "$negative_control_A" ] || [ ! -s "$negative_control_B" ]; then 

        countA=$(( $(wc -l < "$initial_data") ))
        countB=$(( $(wc -l < "$initial_data") ))

        tail -n +2 "$initial_data"  | while IFS=, read -r chr start_gene end_gene lenA lenB strand; do  

            len_region=$(( lenA + lenB ))
            # Temporary files
            sampling_region="$file_name"-sampling-region
            complement_sampling_region="$file_name"-complement-sampling-region

            # Downstream 
            end_sampling_region=$(( end_gene + 5000000 ))
            echo -e "$chr\t$end_gene\t$end_sampling_region" > "$sampling_region"

            intersect "$sampling_region" "$len_region" "$complement_sampling_region"
           
            multiple_exons_negative_control 'downstream' "$complement_sampling_region"
            
            # Upstream 
            start_sampling_region=$(( start_gene - 5000000 ))
            if [ "$start_sampling_region" -lt 0 ]; then start_sampling_region=0; fi 

            echo -e "$chr\t$start_sampling_region\t$start_gene" > "$sampling_region"

            intersect "$sampling_region" "$len_region" "$complement_sampling_region"
           
            multiple_exons_negative_control 'upstream' "$complement_sampling_region" 

			rm -rf "$complement_sampling_region" "$sampling_region" "$random_output" 

        done 

    fi

	## Extract sequences from genome ##
	bedtools getfasta -fi "$genome_seq" -bed "$A_info" -fo "$A_tmp" -s -name -tab  
	bedtools getfasta -fi "$genome_seq" -bed "$B_info" -fo "$B_tmp" -s -name -tab 

	## Reformat file ## 
	reformat_file "$A_tmp" "$A_unsorted"
	reformat_file "$B_tmp" "$B_unsorted"

	sort_output "$A_unsorted" "$negative_control_A" 
	sort_output "$B_unsorted" "$negative_control_B" 

	rm -rf "$A_info" "$B_info" "$A_tmp" "$B_tmp" "$A_unsorted" "$B_unsorted" 
	 
fi

