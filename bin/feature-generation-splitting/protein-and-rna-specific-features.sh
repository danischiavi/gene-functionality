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



############################################################################################################################
initial_data=$1
initial_fasta=$2

first_rna_id=$(awk -F',' 'NR==2 {print $1}' "$initial_data" | tr -d 'RNA')
last_rna_id=$(awk -F',' 'END {print $1}' "$initial_data" | tr -d 'RNA') 

#### File declaration
mkdir -p data/specific && output_directory=data/specific

output_file_interaction="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')interaction.csv"                  # Define name and directory for output file
output_file_MFE="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')MFE.csv" 
output_file_accesibility="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')accesibility.csv" 
output_file_fickett="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')fickett.csv" 
output_file_rnacode="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')RNAcode.csv" 
output_file_rscape="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')covariance.csv" 
output_file_coding_potential="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')coding-potential.csv" 
output_file_structure="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')structure.csv" 

echo InteractionMIN_IntaRNA,InteractionAVE_IntaRNA,InteractionMIN_RNAup,InteractionAVE_RNAup > "$output_file_interaction"
echo MFE > "$output_file_MFE"
echo Accessibility > "$output_file_accesibility"
echo Ficket_score > "$output_file_fickett"
echo RNAcode_score,RNAalifold_score > "$output_file_rnacode"
echo Max_covariance,Min_covariance_Eval > "$output_file_rscape"

#### Create folder for generated r-scape output
rm -rf "$output_directory"/rscapedata && mkdir -p "$output_directory"/rscapedata

######## If multiz100way files already downloaded, use them instead of re-downloading 
if [ -f maf.zip ] 
then 
    unzip data/maf.zip &> /dev/null && rm data/maf.zip
else 
    mkdir -p data/maf  
fi


#### Variables
IntaRNA_exe=IntaRNA
interaction_database=data/raw/curated-interaction-database.fa

RNAfold_exe=RNAfold

access_file=bin/access_py.py
cpc2_file=bin/CPC2_standalone-1.0.1/bin/CPC2.py

############################################################################################################################

# RNA:RNA interactions - Minimum and Average interaction energy 

############################################################################################################################

######## Need to clear any previously set lib path, as otherwise the defined lib path will be appended onto the previous
[ -z "$lib_directory" ] || unset LD_LIBRARY_PATH

var=$first_rna_id                                                                                             
last_seq=$last_rna_id               

# Temporary files
intaRNA_input="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')intaRNA-input"
intaRNA_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')intaRNA-output"
kcal_mol="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')kcal-mol"
RNAup_input="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')RNAup-input"
RNAup_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')RNAup-output"
kcal_mol_RNAup="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')kcal-mol-RNAup"

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
RNAfold_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')RNAfold-output"

#### MFE calculation
$RNAfold_exe < $initial_fasta >> "$RNAfold_output" 2>>errors.log

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
cpc2_output="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')cpc2-output"

## (python script so cannot be added to $PATH)
python3 "$cpc2_file" -i "$initial_fasta" -o "$cpc2_output" >/dev/null 2>>errors.log        

# Extract data
awk 'NR > 1 {printf "%.2f\n", $4}' "$cpc2_output" >> "$output_file_fickett"                             # Ficket score in field 4 of output file

rm -rf "$cpc2_output"

############################################################################################################################

 # RNAalifold scores, MFE, RNAcode score

###########################################################################################################################

#### Directories and temporary files

mkdir -p "$output_directory"/rscape/$(basename "${initial_data%.*}" | sed 's/dataset-//')

maf_directory="$output_directory"/maf/$(basename "${initial_data%.*}" | sed 's/dataset-//')
mkdir -p "$maf_directory"

mafFetch_input="$maf_directory"/mafFetch-input

######## Reformats chromosome coordinates for mafFetch
mafFetch_coordinates="${output_directory}/$(basename "${initial_data%.*}" | sed 's/dataset//')coordinates"
awk -F',' 'NR > 1 {print $3,$4,$5}' "$initial_data" > "$coordinates"


var="$first_rna_id"
#rm -rf *.stk

while IFS=$'\t' read -r line; do

    rna_id=RNA"$var"
    maf_file="$maf_directory"/"$rna_id".maf

    # Obtain MSA files from 100way alignment
    if [ ! -f data/maf/$name/"$rna_id".maf ]; then
        
        echo "$line" > "$mafFetch_input"
        #echo "$mafFetch_exe hg38 multiz100way data/maf/mafFetch-input mafOut 2>>errors.log ;" >> errors.log
	    $mafFetch_exe hg38 multiz100way "$mafFetch_input" "$maf_file" 2>>errors.log ;
	    #do sleep 4 ; done
    fi

    if [ ! -z $maf_file ]; then
        
        ######## Run RNAcode
        RNAcode_output=data/rnacode-output
	    #echo "$RNAcode_exe $file -o $RNAcode_output 2>>errors.log &> /dev/null" >>errors.log
	    $RNAcode_exe $maf_file -o $RNAcode_output 2>>errors.log &> /dev/null
            
        # Takes the largest score, but records zero if no significant hits found.
        rnacode=$( grep -v "No significant coding regions found.\|#\|=" $RNAcode_output | grep . | tr -s ' ' | cut -d ' ' -f 10 | grep . | sort -V | tail -1 )
        [ -z $rnacode ] && rnacode=NA
    
        ####### Convert MSA.maf files to STK format for RNAalifold and R-scape
        maf_to_stockholm  $maf_file $rna_id      # is it better to declare function have it on a different script??  
        STK_input=data/STK/$name/"$rna_id".stk

        ######## Generate secondary structure consensus sequence and associated MFE value
        RNAalifold_output=data/RNAalifold/"$rna_id"-RNAalifold

        echo "$RNAalifold_exe -f S --noPS --aln $STK_input >$rna_id.rnaalifold 2>>errors.log" >>errors.log
        #timeout 60m
	    $RNAalifold_exe -f S --noPS --aln $STK_input > "$RNAalifold_output" 2>>errors.log
        rna_score=$( cat $RNAalifold_output | tail -1 | cut -d ' ' -f 2- | tr -d "(" | cut -d "=" -f 1 )

	    # Remove gap-only columns:
	    #echo "esl-alimask -g --gapthresh 1 data/STK/$rna_id.stk > maf/$rna_id.stk" >>errors.log
	    #esl-alimask -g --gapthresh 1 $rna_id.stk > maf/$rna_id.stk
	    
        [ -z "$rna_score" ] && rna_score=NA
        echo $rnacode,$rna_score | tr -d ' ' >> data/rnacode-$name.csv
	    
        exit_status=$?
	    
        ######## Run R-scape only if previous analyses didn't time out
        Rscape_output=data/rscape/"$rna_id".cov
        if [ "$exit_status" -ne "124" ]; then
        
		    echo "$rscape_exe -E 100 -s $STK_input >/dev/null 2>>errors.log" >>errors.log
            $rscape_exe -E 100 $STK_input > "$Rscape_output" /dev/null 2>>errors.log # deleted -s option since structure is required 
       
        else

            touch $Rscape_output

        fi
    
    else
        ######## Generate blanks if no MAF available
        echo NA,NA >> data/rnacode-$name.csv
        touch $Rscape_output

    fi

    (( var++ ))
    
done < "$mafFetch_coordinates"

zip -r maf maf &> /dev/null
rm -rf maf/
############################################################################################################################

# Process R-scape results to obtain covariance features and associated E-values

############################################################################################################################

######## Process each rscape file individually and in order
count_file="$first_rna_id"
total_rna="$last_rna_id"

while [ $count_file -le $total_rna ]
do
        if [ -f rscapedata/RNA$count_file.cov ] && [ -s rscapedata/RNA$count_file.cov ]  # Check file exists and is not empty
        then
            file=rscapedata/RNA$count_file.cov
            max=$(  grep -r 'GTp' $file       | cut -d '[' -f 2 | tr ']' ',' | cut -d ',' -f 2)  # Max covariance observed
            min=$(  grep -v "#" $file      | cut -f 4 | sort -n | head -1 | tr -d '-' )  # Min covariance observed
            eval=$( grep -m 1 "$min" $file | cut -f 5 )  # E-val of min covariance
            [[ "$eval" == "no significant pairs" ]] && eval=NA   
            echo $max,$eval >> rscape-$name.csv
        else
            ######## If no covariance calculated/no data available
            echo NA,NA >> rscape-$name.csv
        fi
        
        (( count_file++ ))

done

#####?#####
rm -rf score
rm -rf *.ps
rm -rf *.stk
rm -rf coordinates
rm -rf overBed
rm -rf final.stk
rm -rf RNA*.fa
rm -rf info
rm -rf mafOut
rm -rf score
rm -rf mafOut.txt
zip -r rscapedata rscapedata &> /dev/null
rm -rf rscapedata/
########