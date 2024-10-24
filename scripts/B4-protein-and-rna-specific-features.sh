#!/bin/bash
#
# Script Name: protein-and-rna-specific-features.sh
#
#
# Description: This script runs all the scripts to calculates the chosen features for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences. 
#
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.

############################################################################################################################
initial_data=$1
initial_fasta=$2
interaction_database=$3
bigBedToBed_exe=$4
#bigBedToBed_exe=/Volumes/archive/userdata/student_users/danielaschiavinato/dani-scratch/gene-functionality/bin/bigBedToBed
cactus_align_url="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.bigMaf"

first_rna_id=$(awk -F',' 'NR==2 {print $1}' "$initial_data" | tr -d 'RNA')
last_rna_id=$(awk -F',' 'END {print $1}' "$initial_data" | tr -d 'RNA') 

#### File declaration ####
base_name=$(basename "${initial_data%.*}" | sed 's/-dataset//')
mkdir -p data/specific && output_directory=data/specific
file_name="${output_directory}/$base_name"

# Final output file
output_file_specific="$file_name"-specific.csv

# Temporary output files
output_file_interaction="$file_name"-interaction.csv  
output_file_MFE="$file_name"-MFE.csv 
output_file_rnacoding="$file_name"-coding.csv 
output_file_Rscape="$file_name"-covariance.csv


############################################################################################################################

# RNA:RNA interactions - Average interaction energy 

############################################################################################################################

if [ ! -e "$output_file_interaction" ]; then

    echo "Interaction_ave" > "$output_file_interaction"
    
    var=$first_rna_id                                                                                             
    last_seq=$last_rna_id               

    # Temporary files
    RNAup_input="$file_name"-RNAup-input
    RNAup_output="$file_name"-RNAup-output
    kcal_mol_RNAup="$file_name"-kcal-mol-RNAup

    while [ "$var" -le "$last_seq" ]; do

        cat <(grep -w -A 1 ">RNA$var" "$initial_fasta") "$interaction_database" > "$RNAup_input"
    
        RNAup --interaction_first --no_output_file -b --noLP -c 'S' < "$RNAup_input" > "$RNAup_output" 2>>errors.log
        
            # −b, −−include_both : Include the probability of unpaired regions in both (b) RNAs.
            # −−interaction_first : Activate interaction mode using first sequence only.
            # ? −−noLP : Produce structures without lonely pairs (helices of length 1).
            #rm -rf *.out  #no_output_files option is not working.. 

        # Extract min and ave energy from output
        grep '(&)' "$RNAup_output" | awk -F'=' '{print $1}' | awk -F' ' '{print $NF}' | tr -d '(' > "$kcal_mol_RNAup" 
        count_RNAup=$( wc -l < "$kcal_mol_RNAup" )
        min_RNAup=$(head -n 1 "$kcal_mol_RNAup")
        sum_RNAup=0

        while read -r number; do
            ######## If the interaction energy is smaller than the recorded min, update the min
            if (( $(echo "$number < $min_RNAup" | bc -l) )); then min_RNAup=$number; fi
            sum_RNAup=$(echo "scale=3; $sum_RNAup+$number" | bc)
        
        done < "$kcal_mol_RNAup"

        ######## If no interactions were recorded, set minimum and average as NA 
        if [[ "$count_RNAup" == 0 ]]; then
            ave_RNAup=NA
            min_RNAup=NA
        else
            ave_RNAup=$(echo "scale=3; $sum_RNAup/$count_RNAup" | bc )
        fi

        echo "$ave_RNAup" >> "$output_file_interaction"

        (( var++ ))
   
    done

fi

rm -rf "$RNAup_input"
rm -rf "$RNAup_output"
rm -rf "$kcal_mol_RNAup"


############################################################################################################################

# RNA structure 

############################################################################################################################

if [ ! -e "$output_file_MFE" ]; then
    
    ## Temporary file
    echo "MFE" > "$output_file_MFE"
    RNAfold_output="$file_name"RNAfold-output

    #### MFE calculation
    RNAfold < "$initial_fasta" >> "$RNAfold_output" 2>>errors.log

    # Extract the data from output file & replace empty lines with zero                                
    grep "(" "$RNAfold_output" | rev | cut -d "(" -f 1 | rev | tr -d ")" | tr -d " " | awk '!NF{$0="0"}1' >> "$output_file_MFE"      

    rm -rf *.ps 
    rm -rf "$RNAfold_output"

fi

############################################################################################################################

 # RNAalifold scores, MFE, RNAcode score  - RNAalifold; R-scape; RNAcode - MSA: 241way cactus alignment (zoonomia)

###########################################################################################################################

if [ ! -e "$output_file_rnacoding" ] || [ ! -e "$output_file_Rscape" ]; then

    #### Directories and temporary files ####
    
	maf_directory="$output_directory"/maf/"$base_name"
    mkdir -p "$maf_directory"

    RNAalifold_directory="$output_directory"/RNAalifold/"$base_name"
    mkdir -p "$RNAalifold_directory"

    RNAcode_directory="$output_directory"/RNAcode/"$base_name"
    mkdir -p "$RNAcode_directory"

    Rscape_directory="$output_directory"/Rscape/"$base_name"
    mkdir -p  "$Rscape_directory"
    
  
    #### Calculate features line by line #### 
    tail -n +2 "$initial_data" | while IFS=, read -r rna_id _ chr start end _; do

        # Temporary files # 
        Rscape_output="${Rscape_directory}/${rna_id}"
        maf_file="${maf_directory}/${rna_id}.maf"
   
        if [ ! -s "$maf_file" ]; then                                                                               

			#### Obtain MSA maf-files from 241-way cactus alignment ####
            "$bigBedToBed_exe" "$cactus_align_url" stdout -chrom="$chr" -start="$start" -end="$end" | cut -f 4 | tr ';' '\n' > "$maf_file"                                                                   # to extract MSA of regions from UCSC server  
                                                                                        
        fi

        if [ -s "$maf_file" ]; then 
        
            #### Reformat MSA.maf files #### 

            stk_dir="${output_directory}/stk/${base_name}"
			mkdir -p "$stk_dir"
            stk_file="$stk_dir"/"$rna_id".stk

            aln_dir="${output_directory}/aln/${base_name}"
            mkdir -p "$aln_dir"
            aln_file="${aln_dir}/${rna_id}.aln"

            if [ ! -s "$stk_file" ]; then
                
                ## maf to stk for RNAalifold and R-scape
                ./scripts/B4.1-maf-to-stk.sh "$maf_file" "$stk_dir"                                                               # Output: stk file stored at stk_dir 

			fi

			if [ -s "$stk_file" ] && ! grep -q "maf file NA" "$stk_file"; then                                    
                ## stk to clustal (.aln) for Rcode (see notes for more details)
                esl-reformat --mingap clustal "$stk_file" > "$aln_file"
                # --mingap remove columns containing all gaps
				
			else
				echo "maf file was removed - no relevant information" >> errors.log 
				rm -rf "$maf_file" "$stk_file"
			fi
			
			#### Run Features #### 

			##### RNAalifold ####
			if [ -s "$stk_file" ]; then

            	RNAalifold_output="$RNAalifold_directory"/"$rna_id"
				rnaalifold_score=""

            	if [ ! -s "$RNAalifold_output" ]; then                          
                
                	## Generate secondary structure consensus sequence and associated MFE value
                	echo "RNAalifold --input-format=S --noPS $stk_file >$RNAalifold_output 2>>errors.log" >> errors.log
	            	timeout 5m RNAalifold --input-format=S --noPS "$stk_file" > "$RNAalifold_output" 2>> errors.log || echo "RNAalifold_timeout" > "$RNAalifold_output"

					if grep -q "RNAalifold_timeout" "$RNAalifold_output"; then echo "NA" > "$RNAalifold_output"; fi
			
				fi

				rnaalifold_score=$( tail -1 "$RNAalifold_output" | cut -d ' ' -f 2- | tr -d "(" | cut -d "=" -f 1 )
				[ -z "$rnaalifold_score" ] && rnaalifold_score='NA'		
				

			else
				rnaalifold_score=NA
			fi

				
            ##### RNAcode ####
			if [ -s "$aln_file" ]; then 
            
				RNAcode_output="$RNAcode_directory"/"$rna_id"
				rnacode=""

				if [ ! -s "$RNAcode_output" ]; then

	            	echo "RNAcode $aln_file -o $RNAcode_output 2>>errors.log &> /dev/null" >>errors.log
	            	timeout 5m RNAcode "$aln_file" -o "$RNAcode_output" 2>> errors.log &> /dev/null || echo "RNAcode_timeout" > "$RNAcode_output"
			
					if grep -q "RNAcode_timeout" "$RNAcode_output"; then echo "NA" > "$RNAcode_output"; fi
					
				fi
			
				# Takes the largest score, but records zero if no significant hits found.
				rnacode=$( grep -v "No significant coding regions found.\|#\|=" "$RNAcode_output" | grep . | tr -s ' ' | cut -d ' ' -f 10 | grep . | sort -V | tail -1 )
				[ -z "$rnacode" ] && rnacode='NA'
				
            fi

			#### Output scores ####
            if [ ! -e "$output_file_rnacoding" ]; then echo "coding_potential,RNAalifold_score" > "$output_file_rnacoding"; fi
            echo "$rnacode,$rnaalifold_score" | tr -d ' ' >> "$output_file_rnacoding"

            ##### R-scape ####
            if [ -s "$stk_file" ]; then

				if [ ! -s "$Rscape_output" ]; then
                
		        	echo "R-scape -E 100 $stk_file >/dev/null 2>>errors.log" >>errors.log
                	timeout 5m R-scape -E 100 --outdir "$Rscape_directory" --nofigures "$stk_file" >/dev/null 2>>errors.log || touch "$Rscape_output"                           # deleted -s option since structure is required ??

                	# remove extra output files 
                	rm -rf "$Rscape_directory"/"$rna_id".power
                	rm -rf "$Rscape_directory"/"$rna_id".sorted.cov
                fi 

			else
                touch "$Rscape_output"
            fi

        else
            echo "$rna_id multiple alignment file is empty" >> errors.log                                                                                      
           	# Generate blanks if no MAF available
            echo NA,NA >> "$output_file_rnacoding"
            touch "$Rscape_output"
		fi
        
    #rm -rf "$maf_file" "$stk_file" "$aln_file"
	
	done

fi
##### Process R-scape output files ####

if [ ! -s "$output_file_Rscape" ]; then
    
    echo "Max_covariance"> "$output_file_Rscape"

    count_file="$first_rna_id"
    total_rna="$last_rna_id"

    while [ "$count_file" -le "$total_rna" ]; do

        file="$Rscape_directory"/RNA"$count_file".cov

        if [ -s "$file" ]; then                                                                         # Check file exists and is not empty, otherwise NA

            max=$(  grep -r 'GTp' "$file"    | cut -d '[' -f 2 | tr ']' ',' | cut -d ',' -f 2)          # Max covariance observed
        
            echo "$max" >> "$output_file_Rscape"
        
        else
        
            echo "NA" >> "$output_file_Rscape"

        fi
        
        (( count_file++ ))

        rm -rf "$file"

    done

fi


## Join output files for better organization into one
if [ ! -s "$output_file_specific" ]; then

    paste -d',' "$output_file_interaction" "$output_file_rnacoding" "$output_file_Rscape" "$output_file_MFE" > "$output_file_specific"

fi


#### Remove excess files #####

#rm -rf "$output_file_rnacoding"
#rm -rf "$output_file_Rscape" 
#rm -rf "$output_file_MFE"

########