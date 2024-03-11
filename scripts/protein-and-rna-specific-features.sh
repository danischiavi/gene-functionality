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
access_file=$4
cpc2_file=$5
bigBedToBed_exe=$6
cactus_align_url=$7

first_rna_id=$(awk -F',' 'NR==2 {print $1}' "$initial_data" | tr -d 'RNA')
last_rna_id=$(awk -F',' 'END {print $1}' "$initial_data" | tr -d 'RNA') 

#### File declaration ####
mkdir -p data/specific && output_directory=data/specific
file_name="${output_directory}/$(basename "${initial_data%.*}" | sed 's/-dataset//')"

# Final output file
output_file_specific="$file_name"-specific.csv

# Temporary output files
output_file_interaction="$file_name"-interaction.csv  
output_file_MFE="$file_name"-MFE.csv 
output_file_accesibility="$file_name"-accesibility.csv 
output_file_fickett="$file_name"-fickett.csv 
output_file_rnacoding="$file_name"-coding.csv 
output_file_Rscape="$file_name"-covariance.csv

#### Variables ####
RNAup_exe=RNAup
RNAfold_exe=RNAfold
RNAalifold_exe=RNAalifold
RNAcode_exe=RNAcode
Rscape_exe=R-scape
esl_reformat_exe=esl-reformat

############################################################################################################################

# RNA:RNA interactions - Minimum and Average interaction energy 

############################################################################################################################

if [ ! -e "$output_file_interaction" ]; then

    echo "Interaction_min","Interaction_ave" > "$output_file_interaction"
    var=$first_rna_id                                                                                             
    last_seq=$last_rna_id               

    # Temporary files
    RNAup_input="$file_name"-RNAup-input
    RNAup_output="$file_name"-RNAup-output
    kcal_mol_RNAup="$file_name"-kcal-mol-RNAup

    while [ $var -le $last_seq ]; do

        cat <(grep -w -A 1 ">RNA$var" "$initial_fasta") "$interaction_database" > "$RNAup_input"
    
        "$RNAup_exe" --interaction_first --no_output_file -b --noLP -c 'S' < "$RNAup_input" > "$RNAup_output" 2>>errors.log
        
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

        echo "$min_RNAup,$ave_RNAup" >> "$output_file_interaction"

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
    "$RNAfold_exe" < $initial_fasta >> "$RNAfold_output" 2>>errors.log

    # Extract the data from output file & replace empty lines with zero                                
    grep "(" "$RNAfold_output" | rev | cut -d "(" -f 1 | rev | tr -d ")" | tr -d " " | awk '!NF{$0="0"}1' >> "$output_file_MFE"      

    rm -rf *.ps 
    rm -rf "$RNAfold_output"

fi


if [ ! -e "$output_file_accesibility" ]; then 

    echo Accessibility > "$output_file_accesibility"

    #### Accessibility calculation
    for sequence in $( grep -v ">" $initial_fasta ); do

        access=$( timeout 60m python3 "$access_file" -s $sequence 2>>errors.log ) 
        exit_status=$?

        if [ -z "$access" ]; then                           # If no value calculated, record NA
                                                                
            echo 'NA' >> "$output_file_accesibility"

        elif [[ "$access" == "THIS happened"* ]]; then      # If error occurred, record NA
                                                     
            echo 'NA' >> "$output_file_accesibility"
        
        elif [[ "$access" == 'inf' ]]; then                 # If error occurred, record NA
                                                     
            echo 'NA' >> "$output_file_accesibility"

        elif [ "$exit_status" -eq "124" ]; then             # If calculation timed out, record NA
                                                          
            echo 'NA' >> "$output_file_accesibility"

        else
            echo "$access" >> "$output_file_accesibility"

        fi

    done

    sed -i 's/nan/NA/g' "$output_file_accesibility"         # make NA readable by R

fi


############################################################################################################################

# Coding Potential - Fickett Score

############################################################################################################################
if [ ! -s "$output_file_fickett" ]; then

    ## Temporary files
    echo "Fickett_score" > "$output_file_fickett"
    cpc2_output="$file_name"cpc2-output

    ## (python script so cannot be added to $PATH)
    python3 "$cpc2_file" -i "$initial_fasta" -o "$cpc2_output" >/dev/null 2>>errors.log        

    # Extract data
    awk 'NR > 1 {printf "%.2f\n", $4}' "$cpc2_output".txt >> "$output_file_fickett"                             # Ficket score in field 4 of output file

    rm -rf "$cpc2_output".txt

fi

############################################################################################################################

 # RNAalifold scores, MFE, RNAcode score  - RNAalifold; R-scape; RNAcode - MSA: 241way cactus alignment (zoonomia)

###########################################################################################################################

if [ ! -e "$output_file_rnacoding" ] || [ ! -e "$output_file_Rscape" ]; then

    #### Directories and temporary files ####
    maf_directory="$output_directory"/maf/$(basename "${initial_data%.*}" | sed 's/-dataset//')
    mkdir -p "$maf_directory"

    RNAalifold_directory="$output_directory"/RNAalifold/$(basename "${initial_data%.*}" | sed 's/-dataset//')
    mkdir -p "$RNAalifold_directory"

    RNAcode_directory="$output_directory"/RNAcode/$(basename "${initial_data%.*}" | sed 's/-dataset//')
    mkdir -p "$RNAcode_directory"

    Rscape_directory="$output_directory"/Rscape/$(basename "${initial_data%.*}" | sed 's/-dataset//')
    mkdir -p  "$Rscape_directory"
    
  
    #### Obtain MSA maf-files from 241-way cactus alignment ####
    tail -n +2 "$initial_data"  | while IFS=, read -r rna_id _ chr start end _; do

        maf_file="$maf_directory"/"$rna_id".maf

        if [ ! -s "$maf_file" ]; then                                                                               

            "$bigBedToBed_exe" "$cactus_align_url" stdout -chrom="$chr" -start="$start" -end="$end" | cut -f 4 | tr ';' '\n' > "$maf_file"                                                                   # to extract MSA of regions from UCSC server  
                                                                                        
        fi

        if [ -s "$maf_file" ]; then 
        
            #### Reformat MSA.maf files for tools #### 

            stk_dir="${output_directory}/stk/$(basename "${initial_data%.*}" | sed 's/-dataset//')"
            if [ ! -e "$stk_dir" ]; then mkdir -p "$stk_dir"; fi
            stk_file="$stk_dir"/"$rna_id".stk

            aln_dir="${output_directory}/aln/$(basename "${initial_data%.*}" | sed 's/-dataset//')"
            mkdir -p "$aln_dir"
            aln_file="$aln_dir"/"$rna_id".aln

            if [ ! -s "$stk_file" ]; then
                
                ## maf to stk for RNAalifold and R-scape
                ./scripts/maf-to-stk.sh "$maf_file" "$stk_dir"                                                               # Output: stk file stored at stk_dir 
                                                    
                ## stk to clustal (.aln) for Rcode (see notes for more details)
                "$esl_reformat_exe" --mingap clustal "$stk_file" > "$aln_file"
                # --mingap remove columns containing all gaps
            
            fi

            ##### RNAalifold ####
            RNAalifold_output="$RNAalifold_directory"/"$rna_id"

            if [ ! -s "$RNAalifold_output" ]; then                          
                
                ## Generate secondary structure consensus sequence and associated MFE value
                echo "$RNAalifold_exe --input-format=S --noPS $stk_file >$RNAalifold_output 2>>errors.log" >>errors.log
                #timeout 60m
	            "$RNAalifold_exe" --input-format=S --noPS "$stk_file" > "$RNAalifold_output" 2>>errors.log
            
            fi

            rnaalifold_score=$( cat "$RNAalifold_output" | tail -1 | cut -d ' ' -f 2- | tr -d "(" | cut -d "=" -f 1 )
            
            ##### RNAcode ####
            RNAcode_output="$RNAcode_directory"/"$rna_id"
            if [ ! -s "$RNAcode_output" ]; then

	            echo "$RNAcode_exe $aln_file -o $RNAcode_output 2>>errors.log &> /dev/null" >>errors.log
	            "$RNAcode_exe" "$aln_file" -o "$RNAcode_output" 2>>errors.log &> /dev/null
            
            fi
            # Takes the largest score, but records zero if no significant hits found.
            rnacode=$( grep -v "No significant coding regions found.\|#\|=" "$RNAcode_output" | grep . | tr -s ' ' | cut -d ' ' -f 10 | grep . | sort -V | tail -1 )
            
            [ -z $rnacode ] && rnacode=NA
            [ -z "$rnaalifold_score" ] && rnaalifold_score=NA
            
            if [ ! -e "$output_file_rnacoding" ]; then echo "RNAcode_score,RNAalifold_score" > "$output_file_rnacoding"; fi
            echo "$rnacode,$rnaalifold_score" | tr -d ' ' >> "$output_file_rnacoding"

            ##### R-scape ####
            Rscape_output="$Rscape_directory"/"$rna_id"
            if [ ! -s "$Rscape_output" ]; then
                
		        echo "$Rscape_exe -E 100 $stk_file >/dev/null 2>>errors.log" >>errors.log
                "$Rscape_exe" -E 100 --outdir "$Rscape_directory" --nofigures "$stk_file" >/dev/null 2>>errors.log                           # deleted -s option since structure is required ??

                # remove extra output files 
                rm -rf "$Rscape_directory"/"$rna_id".power
                rm -rf "$Rscape_directory"/"$rna_id".sorted.cov
                 
            fi

        else                                                                                      # Generate blanks if no MAF available
           
            echo NA,NA >> "$output_file_rnacoding"
            touch "$Rscape_output"

        fi

    done 

fi
##### Process R-scape output files ####

if [ ! -s "$output_file_Rscape" ]; then
    
    echo Max_covariance,Min_covariance_Eval > "$output_file_Rscape"

    count_file="$first_rna_id"
    total_rna="$last_rna_id"

    while [ $count_file -le $total_rna ]; do

        file="$Rscape_directory"/RNA"$count_file".cov

        if [ -s "$file" ]; then                                                                         # Check file exists and is not empty, otherwise NA

            max=$(  grep -r 'GTp' "$file"    | cut -d '[' -f 2 | tr ']' ',' | cut -d ',' -f 2)          # Max covariance observed
            min=$(  grep -v "#" "$file" | grep -v '^$' | cut -f 4 | sort -n | head -n 1 | tr -d '-' )   # Min covariance observed
            eval=$( grep -m 1 "$min" "$file" | cut -f 5 )                                               # E-val of min covariance
            
            [[ "$eval" == "no significant pairs" ]] && eval=NA   
            
            echo "$max,$eval" >> "$output_file_Rscape"
        
        else
        
            echo NA,NA >> "$output_file_Rscape"

        fi
        
        (( count_file++ ))

    done

fi


## Join output files for better organization into one
if [ ! -s "$output_file_specific" ]; then
    paste -d',' "$output_file_interaction" "$output_file_fickett" "$output_file_rnacoding" "$output_file_accesibility" "$output_file_Rscape" "$output_file_MFE" > "$output_file_specific"
fi


#### Remove excess files #####

#rm -rf "$output_file_fickett"
#rm -rf "$output_file_rnacoding"
#rm -rf "$output_file_accesibility" 
#rm -rf "$output_file_Rscape" 
#rm -rf "$output_file_MFE"

########