#!/bin/bash
#
# Script Name: protein-and-rna-specific-features.sh
#
#
# Description: This script runs all the scripts to calculates the chosen features for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences. 
#
# 
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.

############################################################################################################################
initial_data=$1
initial_fasta=$2

first_rna_id=$(awk -F',' 'NR==2 {print $1}' "$initial_data" | tr -d 'RNA')
last_rna_id=$(awk -F',' 'END {print $1}' "$initial_data" | tr -d 'RNA') 

#### File declaration ####
mkdir -p data/specific && output_directory=data/specific
file_name="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')"

output_file_interaction="$file_name"interaction.csv                  # Define name and directory for output file
output_file_MFE="$file_name"MFE.csv 
output_file_accesibility="$file_name"accesibility.csv 
output_file_fickett="$file_name"fickett.csv 
output_file_rnacoding="$file_name"RNAcoding.csv 
output_file_rscape="$file_name"covariance.csv

output_file_coding_potential="$file_name"coding-potential.csv 
output_file_structure="$file_name"structure.csv 

echo InteractionMIN_IntaRNA,InteractionAVE_IntaRNA,InteractionMIN_RNAup,InteractionAVE_RNAup > "$output_file_interaction"
echo MFE > "$output_file_MFE"
echo Accessibility > "$output_file_accesibility"
echo Ficket_score > "$output_file_fickett"
echo RNAcode_score,RNAalifold_score > "$output_file_rnacoding"
echo Max_covariance,Min_covariance_Eval > "$output_file_rscape"


#### Variables ####
IntaRNA_exe=IntaRNA
interaction_database=data/raw/curated-interaction-database.fa

RNAfold_exe=RNAfold

access_file=bin/access_py.py
cpc2_file=bin/CPC2_standalone-1.0.1/bin/CPC2.py

bigBedToBed_exe=bin/bigBedToBed
cactus_align_url=http://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.bigMaf

RNAcode_exe=RNAcode
Rscape_exe=R-scape
esl_reformat_exe=esl-reformat

############################################################################################################################

# RNA:RNA interactions - Minimum and Average interaction energy 

############################################################################################################################

######## Need to clear any previously set lib path, as otherwise the defined lib path will be appended onto the previous
[ -z "$lib_directory" ] || unset LD_LIBRARY_PATH

var=$first_rna_id                                                                                             
last_seq=$last_rna_id               

# Temporary files
intaRNA_input="$file_name"intaRNA-input
intaRNA_output="$file_name"intaRNA-output
kcal_mol="$file_name"kcal-mol
RNAup_input="$file_name"RNAup-input
RNAup_output="$file_name"RNAup-output
kcal_mol_RNAup="$file_name"kcal-mol-RNAup

while [ $var -le $last_seq ]; do

    echo "$( grep -w -A 1 ">RNA$var" $initial_fasta )" > "$intaRNA_input"
    
    #if [ -z "$lib_directory" ]                                                                               # If Boost library didn't have to specified, run as normal.
    #then
    $IntaRNA_exe -t "$interaction_database" -m "$intaRNA_input" > "$intaRNA_output" 2>>errors.log           
    #else
    #    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$lib_directory $IntaRNA_exe -m $intaRNA_input -t $interaction_database > data/intaRNA-results 2>>errors.log
    #fi 
    
    ######## Grab all recorded interactions
  
    grep "energy" "$intaRNA_output" | cut -d ':' -f 2 | tr -d "kcal/mol" | tr -d ' ' > "$kcal_mol"        
                                                                                                      
    count=$(wc -l < "$kcal_mol" )                                                                        # Number of interaction energies recorded
    min=$( head -n 1 "$kcal_mol")                                                                        # Min interaction energy
    sum=0                                                                                                # Sum of all energies calculated, prior to averaging

    while read -r number; do 
        
        if (( $(echo "$number < $min" | bc -l) )); then min="$number"; fi                                   # If the interaction energy is smaller than the recorded min, update the min
        sum=$(echo "scale=3; $sum+$number" | bc)
    
    done < "$kcal_mol"

    ######## If no interactions were recorded, set minimum and average as NA 
    if [[ "$count" == 0 ]]; then
        ave='NA'
        min='NA'
    else
        ave=$(echo "scale=3; $sum/$count" | bc )
    fi

    #### RNAup

    # RNA input: sequence following by the curated sequences
    cat "$intaRNA_input" "$interaction_database" > "$RNAup_input"

    # Run
    RNAup --interaction_first --no_output_file -b --noLP -c 'S' < "$RNAup_input" > "$RNAup_output"
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

    echo "$min,$ave,$min_RNAup,$ave_RNAup" >> "$output_file_interaction"

    (( var++ ))
   
done

rm -rf "$intaRNA_input"
rm -rf "$intaRNA_output"
rm -rf "$kcal_mol"
rm -rf "$RNAup_input"
rm -rf data/intaRNA-output
rm -rf "$RNAup_output"
rm -rf "$kcal_mol_RNAup"


############################################################################################################################

# RNA structure 

############################################################################################################################

## Temporary file
RNAfold_output="$file_name"RNAfold-output

#### MFE calculation
"$RNAfold_exe" < $initial_fasta >> "$RNAfold_output" 2>>errors.log

# Extract the data from output file & replace empty lines with zero                                
grep "(" "$RNAfold_output" | rev | cut -d "(" -f 1 | rev | tr -d ")" | tr -d " " | awk '!NF{$0="0"}1' >> "$output_file_MFE"      

rm -rf *.ps 
rm -rf "$RNAfold_output"

   
#### Accessibility calculation
for sequence in $( grep -v ">" $initial_fasta ); do

    access=$( timeout 60m python3 "$access_file" -s $sequence 2>>errors.log ) 
    exit_status=$?

    if [ -z "$access" ]; then                           # If no value calculated, record NA
                                                                
        echo 'NA' >> "$output_file_accesibility"

    elif [[ "$access" == "THIS happened"* ]]; then      # If error occurred, record NA
                                                     
        echo 'NA' >> "$output_file_accesibility"

    elif [ "$exit_status" -eq "124" ]; then             # If calculation timed out, record NA
                                                          
        echo 'NA' >> "$output_file_accesibility"

    else
        echo "$access" >> "$output_file_accesibility"

    fi

done

sed -i 's/nan/NA/g' "$output_file_accesibility"         # make NA readable by R


############################################################################################################################

# Coding Potential - Ficket Score

############################################################################################################################

## Temporary files
cpc2_output="$file_name"cpc2-output

## (python script so cannot be added to $PATH)
python3 "$cpc2_file" -i "$initial_fasta" -o "$cpc2_output" >/dev/null 2>>errors.log        

# Extract data
awk 'NR > 1 {printf "%.2f\n", $4}' "$cpc2_output" >> "$output_file_fickett"                             # Ficket score in field 4 of output file

rm -rf "$cpc2_output"

############################################################################################################################

 # RNAalifold scores, MFE, RNAcode score  - RNAalifold; R-scape; RNAcode - MSA: 241way cactus alignment (zoonomia)

###########################################################################################################################

#### Directories and temporary files

Rscape_directory="$output_directory"/rscape/$(basename "${initial_data%.*}" | sed 's/dataset-//')
mkdir -p  "$rscape_directory"

maf_directory="$output_directory"/maf/$(basename "${initial_data%.*}" | sed 's/dataset-//')
mkdir -p "$maf_directory"

RNAcode_output="$file_name"rnacode-output
RNAalifold_output="${output_directory}/RNAalifold/$(basename "${initial_data%.*}" | sed 's/dataset//')/$rna_id"
Rscape_output="$Rscape_directory"/"$rna_id".cov


#### Obtain MSA maf-files from 241-way cactus alignment
{
    read  
    while IFS=, read -r rna_id _ chr start end _; do

        maf_file="$maf_directory"/"$rna_id".maf

        if [ ! -e "$maf_file" ]; then                                                                               

            "$bigBedToBed_exe" "$cactus_align_url" \                                                                # to extract MSA of regions from UCSC server  
            stdout -chrom="$chr" -start="$start" -end="$end" | \
            cut -f 4 | tr ';' \                                                                                     # format bed output to maf
            '\n' > "$maf_file"

        fi

        if [ ! -z "$maf_file" ]; then
        
            #### Reformat MSA.maf files for tools #### 

                ## maf to stk for RNAalifold and R-scape

            stk_dir="${output_directory}/stk/$(basename "${initial_data%.*}" | sed 's/-dataset//')"
            mkdir -p "$stk_dir"
            stk_file="$stk_dir"/"$rna_id".stk

            ./bin/maf_to_stk.sh "$maf_file" "$stk_dir"                                                               # Output: stk file stored at stk_dir 
                                                         
            
                ## stk to clustal (.aln) for Rcode (see notes for more details)
            
            aln_dir="${output_directory}/aln/$(basename "${initial_data%.*}" | sed 's/-dataset//')"
            mkdir -p "$aln_dir"
            aln_file="$aln_dir"/"$rna_id".aln

            "$esl_reformat_exe" clustal "$stk_file" > "$aln_file"


            ##### RNAalifold ####

            ## Generate secondary structure consensus sequence and associated MFE value
            
            echo "$RNAalifold_exe -f S --noPS --aln $stk_file >$RNAalifold_output 2>>errors.log" >>errors.log
            #timeout 60m
	        "$RNAalifold_exe" -f S --noPS --aln "$stk_file" > "$RNAalifold_output" 2>>errors.log
            rnaalifold_score=$( cat "$RNAalifold_output" | tail -1 | cut -d ' ' -f 2- | tr -d "(" | cut -d "=" -f 1 )

	        # Remove gap-only columns:
	        #echo "esl-alimask -g --gapthresh 1 data/STK/$rna_id.stk > maf/$rna_id.stk" >>errors.log
	        #esl-alimask -g --gapthresh 1 $rna_id.stk > maf/$rna_id.stk
            exit_status=$?


            ##### R-scape ####

            ## Run only if previous analyses didn't time out
            
            if [ "$exit_status" -ne "124" ]; then
        
		        echo "$Rscape_exe -E 100 -s $stk_file >/dev/null 2>>errors.log" >>errors.log
                "$Rscape_exe" -E 100 "$stk_file" > "$Rscape_output" /dev/null 2>>errors.log # deleted -s option since structure is required 
       
            else

                touch "$Rscape_output"

            fi


            ##### RNAcode ####
            
	        echo "$RNAcode_exe $aln_file -o $RNAcode_output 2>>errors.log &> /dev/null" >>errors.log
	        "$RNAcode_exe" "$aln_file" -o "$RNAcode_output" 2>>errors.log &> /dev/null
            
            # Takes the largest score, but records zero if no significant hits found.
            rnacode=$( grep -v "No significant coding regions found.\|#\|=" "$RNAcode_output" | grep . | tr -s ' ' | cut -d ' ' -f 10 | grep . | sort -V | tail -1 )
            
            ###################

            ## Scores to output file (R-scape output files process below)

            [ -z $rnacode ] && rnacode=NA
            [ -z "$rnaalifold_score" ] && rnaalifold_score=NA

            echo "$rnacode,$rnaalifold_score" | tr -d ' ' >> "$output_file_rnacodig"

        else
            ## Generate blanks if no MAF available
            echo NA,NA >> "$output_file_rnacoding"
            touch "$Rscape_output"

        fi
    
    done 

} < "$initial_data"


##### Process R-scape output files ####

count_file="$first_rna_id"
total_rna="$last_rna_id"

while [ $count_file -lt $total_rna ]; do

    file="$Rscape_directory"/"$count_file".cov

    if [ -f "$file" ] && [ -s "$file" ]; then                                                         # Check file exists and is not empty, otherwise NA

        max=$(  grep -r 'GTp' "$file"    | cut -d '[' -f 2 | tr ']' ',' | cut -d ',' -f 2)            # Max covariance observed
        min=$(  grep -v "#" "$file"      | cut -f 4 | sort -n | head -1 | tr -d '-' )                 # Min covariance observed
        eval=$( grep -m 1 "$min" "$file" | cut -f 5 )                                                 # E-val of min covariance
            
        [[ "$eval" == "no significant pairs" ]] && eval=NA   
            
        echo "$max,$eval" >> "$output_file_rscape"
        
    else
        
        echo NA,NA >> "$output_file_rscape"

    fi
        
    (( count_file++ ))

done

## Join some output files for better organization 
paste -d',' "$output_file_fickett" "$output_file_rnacoding" > "$output_file_coding_potential"
paste -d',' "$output_file_accesibility" "$output_file_rscape" "$output_file_MFE"  > "$output_file_structure"


#####?#####
#rm -rf score
#rm -rf *.ps
#rm -rf *.stk
#rm -rf overBed

#rm -rf RNA*.fa
#rm -rf info
#rm -rf mafOut
#rm -rf score
#rm -rf mafOut.txt
#zip -r rscapedata rscapedata &> /dev/null
#rm -rf rscapedata/
#zip -r maf maf &> /dev/null
#rm -rf maf/
########