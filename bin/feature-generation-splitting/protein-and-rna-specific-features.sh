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

####Checks 

## Interaction database for IntaRNA

if [ -f $additional_folder/curated-interaction-database.fa ] 
then
    interaction_database=$additional_folder/curated-interaction-database.fa
else
    until [ -f $interaction_database ] ; do read -p "Please enter custom name of RNA:RNA interaction database (fa): " i ; interaction_database=$additional_folder/$i ; echo ; done
fi

## IntaRNA

if find $additional_folder/ -executable -type f | grep -q IntaRNA
then
    IntaRNA_exe=$additional_folder/IntaRNA
elif command -v IntaRNA &> /dev/null
then
    IntaRNA_exe=IntaRNA
else
    echo Please check that IntaRNA is installed.
    exit 1
fi

intaRNA_error=$( $IntaRNA_exe -h 2>&1 )  # Need to find out if Boost library needs to be defined separately for IntaRNA to run

if [ ! -z "$error" ] && [[ "$error" != "IntaRNA predictors RNA-RNA interactions"* ]]
then
    read -p "Enter path for Boost C++ lib folder for IntaRNA: " lib_directory
    echo
    if [ -d $lib_directory/ ] ; then : ; else echo Please check that the Boost C++ lib directory is correct. ; exit 1 ; fi
else
    :
fi

## mafFetch

if find $additional_folder/ -executable -type f | grep -q mafFetch
then
    mafFetch_exe=$additional_folder/mafFetch
elif command -v mafFetch &> /dev/null
then
    mafFetch_exe=mafFetch
else
    echo Please check that mafFetch is installed.
    exit 1
fi

## RNAcode

if find $additional_folder/ -executable -type f | grep -q RNAcode
then
    RNAcode_exe=$additional_folder/RNAcode
elif command -v RNAcode &> /dev/null
then
    RNAcode_exe=RNAcode
else
    echo Please check that RNAcode is installed.
    exit 1
fi

## RNAalifold

if find $additional_folder/ -executable -type f | grep -q RNAalifold
then
    RNAalifold_exe=$additional_folder/RNAalifold
elif command -v RNAalifold &> /dev/null
then
    RNAalifold_exe=RNAalifold
else
    echo Please check that RNAalifold is installed.
    exit 1
fi

## RNAfold

if find $additional_folder/ -executable -type f | grep -q RNAfold
then
    RNAfold_exe=$additional_folder/RNAfold
elif command -v RNAfold &> /dev/null
then
    RNAfold_exe=RNAfold
else
    echo Please check that RNAfold is installed.
    exit 1
fi

## RNAplfold

if command -v RNAplfold &> /dev/null
then
    :
else
    echo Please check that RNAplfold is installed and exported to $PATH.
    exit 1
fi

## R-scape

if find $additional_folder/ -executable -type f | grep -q 'R-scape'
then
    rscape_exe=$additional_folder/R-scape
elif command -v R-scape &> /dev/null
then
    rscape_exe=R-scape
else
    echo Please check that R-scape is installed.
    exit 1
fi

## access_py.py

if [ -f $additional_folder/access_py.py ]
then
    :
else
    echo Please check that access_py.py is available in $additional_folder.
fi

## CPC2

if [ -f $additional_folder/CPC2.py ]
then
    cpc2_directory=$additional_folder/CPC2.py
else
    read -p "Enter path to installed CPC2 bin: " cpc2_directory
    echo
    if [ -d $cpc2_directory/ ] ; then : ; else echo Please check that CPC2 is installed. ; exit 1 ; fi
fi


############################################################################################################################

#### File declaration
echo InteractionMIN_IntaRNA,InteractionAVE_IntaRNA,InteractionMIN_RNAup,InteractionAVE_RNAup > data/interaction-intermediate-$name.csv
echo MFE > data/MFE-$name.csv
echo Accessibility > data/access-$name.csv
echo Ficket_score > data/CPC2-"$name".csv
echo RNAcode_score,RNAalifold_score > data/rnacode-$name.csv
echo Max_covariance,Min_covariance_Eval > data/rscape-$name.csv

#### Create folder for generated r-scape output
rm -rf data/rscapedata && mkdir -p data/rscapedata

######## If multiz100way files already downloaded, use them instead of re-downloading 
if [ -f maf.zip ] 
then 
    unzip data/maf.zip &> /dev/null && rm data/maf.zip
else 
    mkdir -p data/maf  
fi

############################################################################################################################

# RNA:RNA interactions - Minimum and Average interaction energy 

############################################################################################################################

######## Need to clear any previously set lib path, as otherwise the defined lib path will be appended onto the previous
[ -z "$lib_directory" ] || unset LD_LIBRARY_PATH

var=$first_rna_id                                                                                             # Global variables for given id number to first and last RNA sequence of the dataset
last_seq=$last_rna_id               

while [ $var -le $last_seq ]
do
    echo "$( grep -w -A 1 ">RNA$var" $initial_fasta )" > data/intaRNA-input
    
    
    #if [ -z "$lib_directory" ]                                                                               # If Boost library didn't have to specified, run as normal.
    #then
    $IntaRNA_exe -t $interaction_database -m data/intaRNA-input > data/intaRNA-output 2>>errors.log           #swap -m for -q (its how IntaRNA works in my session...)
    #else
    #    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$lib_directory $IntaRNA_exe -m $intaRNA_input -t $interaction_database > data/intaRNA-results 2>>errors.log
    #fi 
    
    ######## Grab all recorded interactions
    grep "energy" data/intaRNA-output | cut -d ':' -f 2 | tr -d "kcal/mol" | tr -d ' ' > data/kcal_mol        #It was: grep "interaction energy" intaRNA-results | cut -d '=' -f 2 | tr -d "kcal/mol" | tr -d ' ' > kcal_mol --> but for me output -> energy: -5.20598 kcal/mol

                                                                                                              # Min interaction energy
    count=$(wc -l < data/kcal_mol )                                                                           # Number of interaction energies recorded
    min=$( head -n 1 data/kcal_mol) 
    sum=0                                                                                                     # Sum of all energies calculated, prior to averaging

    while read -r number; do 
        ######## If the interaction energy is smaller than the recorded min, update the min
        if (( $(echo "$number < $min" | bc -l) )); then
            min=$number 
        fi

        sum=$(echo "scale=3; $sum+$number" | bc)
    
    done < data/kcal_mol

    ######## If no interactions were recorded, set minimum and average as NA 
    if [[ "$count" == 0 ]]; then
        ave='NA'
        min='NA'
    else
        ave=$(echo "scale=3; $sum/$count" | bc )
    fi

    #### RNAup

    # RNA input: sequence following by the curated sequences
    cat data/intaRNA-input data/raw/curated-interaction-database.fa > data/RNAup-input

    # Run
    RNAup --interaction_first --no_output_file -b --noLP -c 'S' < data/RNAup-input > data/RNAup.out
    # −b, −−include_both : Include the probability of unpaired regions in both (b) RNAs.
    # −−interaction_first : Activate interaction mode using first sequence only.
    # ? −−noLP : Produce structures without lonely pairs (helices of length 1).
    rm -rf *.out  #no_output_files option is not working.. 

    # Extract min and ave energy from output
    grep '(&)' data/RNAup.out | awk -F'=' '{print $1}' | awk -F' ' '{print $NF}' | tr -d '(' > data/kcal_mol-RNAup 
    count_RNAup=$( wc -l < data/kcal_mol-RNAup )
    min_RNAup=$(head -n 1 data/kcal_mol-RNAup)
    sum_RNAup=0

    while read -r number; do
        ######## If the interaction energy is smaller than the recorded min, update the min
        if (( $(echo "$number < $min_RNAup" | bc -l) )); then
            min_RNAup=$number
        fi
        
        sum_RNAup=$(echo "scale=3; $sum_RNAup+$number" | bc)
        
    done < data/kcal_mol-RNAup

    ######## If no interactions were recorded, set minimum and average as NA 
    if [[ "$count_RNAup" == 0 ]]; then
        ave_RNAup=NA
        min_RNAup=NA
    else
        ave_RNAup=$(echo "scale=3; $sum_RNAup/$count_RNAup" | bc )
    fi


    echo RNA$var,$min,$ave,$min_RNAup,$ave_RNAup >> data/interaction-intermediate-$name.csv

    (( var++ ))
   
done

rm -rf data/intaRNA-results
rm -rf data/kcal_mol
rm -rf data/kcal_mol-RNAup
rm -rf data/intaRNA-input
rm -rf data/intaRNA-output
rm -rf data/RNAup-input
rm -rf data/RNAup-output
rm -rf data/RNAup.out

############################################################################################################################

# RNA structure 

############################################################################################################################

#### MFE calculation
$RNAfold_exe < $initial_fasta >> data/rnafold-output 2>>errors.log

# Extract the data from output file
grep "(" data/rnafold-output | rev | cut -d "(" -f 1 | rev | tr -d ")" | tr -d " " >> data/MFE.csv
awk '!NF{$0="0"}1' data/MFE.csv >> data/MFE-$name.csv                                               # Replace empty lines with zero
rm -rf *.ps
rm -rf data/MFE.csv 
rm -rf data/rnafold-output


#### Accessibility calculation
for sequence in $( grep -v ">" $initial_fasta )
do 
    access=$( timeout 60m python3 bin/access_py.py -s $sequence 2>>errors.log ) 
    exit_status=$?
    if [ -z "$access" ]; then                                                                       # If no value calculated, record NA
        echo 'NA' >> data/access-$name.csv
    elif [[ "$access" == "THIS happened"* ]]; then                                                  # If error occurred, record NA
        echo 'NA' >> data/access-$name.csv
    elif [ "$exit_status" -eq "124" ]; then                                                         # If calculation timed out, record NA
        echo 'NA' >> data/access-$name.csv
    else
        echo "$access" >> data/access-$name.csv
    fi
done

sed -i 's/nan/NA/g' data/access-$name.csv                                                           # make NA readable by R


############################################################################################################################

# Coding Potential - Ficket Score

############################################################################################################################

## (python script so cannot be added to $PATH)

python3 bin/CPC2_standalone-1.0.1/bin/CPC2.py -i $initial_fasta -o data/cpc2-output >/dev/null 2>>errors.log        

# Extract data
awk 'NR > 1 {printf "%.2f\n", $4}' data/cpc2-output.txt >> data/CPC2-"$name".csv                                    # Ficket score in field 4 of output file, extract value with 2 decimals
rm -rf data/cpc2-output.txt

############################################################################################################################

 # RNAalifold scores, MFE, RNAcode score

###########################################################################################################################

# Obtain MSA files from 241-way cactus alignment

bigBedToBed_exe='bin/bigBedToBed'       # to extract MSA of regions from UCSC server  
241way_url='http://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.bigMaf'
input_file='data/initial_data_sample'

{
    read  
    while IFS=, read -r rna_id _ chr start end _; do

        # Obtain MSA files from 241-way cactus alignment
        maf_file='data/zoonomia-maf/$rna_id.maf'
        $bigBedToBed_exe $241way_url \
        stdout -chrom=$chr -start=$start -end=$end | cut -f 4 | tr ';' \        # format bed output to maf
        '\n' > $maf_file

    if [ ! -z $maf_file ]; then
        
        ######## Run RNAcode
        RNAcode_output='data/rnacode-output'
	    #echo "$RNAcode_exe $file -o $RNAcode_output 2>>errors.log &> /dev/null" >>errors.log
	    $RNAcode_exe $maf_file -o $RNAcode_output 2>>errors.log &> /dev/null
            
        # Takes the largest score, but records zero if no significant hits found.
        rnacode=$( grep -v "No significant coding regions found.\|#\|=" $RNAcode_output | grep . | tr -s ' ' | cut -d ' ' -f 10 | grep . | sort -V | tail -1 )
        [ -z $rnacode ] && rnacode=NA
    
        ####### Convert MSA.maf files to STK format for RNAalifold and R-scape
        bin/maf_to_stockholm  $maf_file $rna_id      # is it better to declare function have it on a different script??  
        STK_input='data/STK/$rna_id.stk'

        ######## Generate secondary structure consensus sequence and associated MFE value
        RNAalifold_output='data/RNAalifold/"$rna_id"alifold'

        echo "$RNAalifold_exe -f S --noPS --aln $STK_input >$rna_id.rnaalifold 2>>errors.log" >>errors.log
        #timeout 60m
	    $RNAalifold_exe -f S --noPS --aln $STK_input > $RNAalifold_output 2>>errors.log
        rna_score=$( cat $RNAalifold_output | tail -1 | cut -d ' ' -f 2- | tr -d "(" | cut -d "=" -f 1 )

	    # Remove gap-only columns:
	    #echo "esl-alimask -g --gapthresh 1 data/STK/$rna_id.stk > maf/$rna_id.stk" >>errors.log
	    #esl-alimask -g --gapthresh 1 $rna_id.stk > maf/$rna_id.stk
	    
        [ -z "$rna_score" ] && rna_score=NA
        echo $rnacode,$rna_score | tr -d ' ' >> data/rnacode-$name.csv
	    
        exit_status=$?
	    
        ######## Run R-scape only if previous analyses didn't time out
        Rscape_output='data/rscape/"$rna_id".cov'
        if [ "$exit_status" -ne "124" ]; then
        
		    echo "$rscape_exe -E 100 -s $STK_input >/dev/null 2>>errors.log" >>errors.log
            $rscape_exe -E 100 -s $STK_input > $Rscape_output /dev/null 2>>errors.log
       
        else

            touch $Rscape_output

        fi
    
    else
        ######## Generate blanks if no MAF available
        echo NA,NA >> data/rnacode-$name.csv
        touch $Rscape_output

    fi
    
    done 

} < $initial_data

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