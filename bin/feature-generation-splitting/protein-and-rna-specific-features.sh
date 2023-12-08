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
echo InteractionMIN,InteractionAVE > ./data/interaction-intermediate.csv
echo MFE > ./data/MFE.csv
echo Accessibility > ./data/access.csv
echo RNAcode_score,RNAalifold_score > ./data/rnacode.csv
echo Max_covariance,Min_covariance_Eval > ./data/rscape-dataset.csv

#### Create folder for generated r-scape output
rm -rf ./data/rscapedata && mkdir -p ./data/rscapedata

######## If multiz100way files already downloaded, use them instead of re-downloading 
if [ -f maf.zip ] 
then 
    unzip ./data/maf.zip &> /dev/null && rm ./data/maf.zip
else 
    mkdir -p ./data/maf  
fi

############################################################################################################################

# RNA:RNA interactions - Minimum and Average interaction energy 

############################################################################################################################

######## Need to clear any previously set lib path, as otherwise the defined lib path will be appended onto the previous
[ -z "$lib_directory" ] || unset LD_LIBRARY_PATH

for seq in $( grep ">" $initial_fasta )
do
    grep -A 1 $seq $initial_fasta >target.fa 
    if [ -z "$lib_directory" ]   # If Boost library didn't have to specified, run as normal.
    then
        $IntaRNA_exe -t $interaction_database -m target.fa > intaRNA-results 2>>errors.log #swap -m for -q (its how IntaRNA works in my session...)
    else
        LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$lib_directory $IntaRNA_exe -q target.fa -t $interaction_database > intaRNA-results 2>>errors.log
    fi 
    
    ######## Grab all recorded interactions
    grep "energy" intaRNA-results | cut -d ':' -f 2 | tr -d "kcal/mol" | tr -d ' ' > kcal_mol #It was: grep "interaction energy" intaRNA-results | cut -d '=' -f 2 | tr -d "kcal/mol" | tr -d ' ' > kcal_mol --> but for me output -> energy: -5.20598 kcal/mol

    min=0    # Min interaction energy
    count=$( grep -c "energy" intaRNA-results )   # Number of interaction energies recorded
    sum=0    # Sum of all energies calculated, prior to averaging

    for number in $( cat kcal_mol )
    do
        ######## If the interaction energy is smaller than the recorded min, update the min
        if (( $( echo "$number < $min" | bc -l ) ))
        then
            min=$number
            sum=$(echo "scale=3; $sum+$number" | bc)
        else
            sum=$(echo "scale=3; $sum+$number" | bc)
        fi
    done

    ######## If no interactions were recorded, set minimum and average as NA 
    if (( $( echo "$count == 0" | bc -l) ))
    then
        ave=NA
        min=NA
    else
        ave=$(echo "scale=3; $sum/$count" | bc )
    fi

    echo $min,$ave >> interaction-intermediate.csv

done

rm -rf target.fa
rm -rf intaRNA-results
rm -rf kcal_mol


############################################################################################################################

# RNA structure 

############################################################################################################################

## MFE calculation
$RNAfold_exe < $initial_fasta >> rnafold-output 2>>errors.log

grep "(" rnafold-output | rev | cut -d "(" -f 1 | rev | tr -d ")" | tr -d " " >> ./data/MFE.csv
awk '!NF{$0="0"}1' ./data/MFE.csv > ./data/MFE-final.csv
rm -rf *.ps
rm -fr ./data/MFE.csv 


## Accessibility calculation
for sequence in $( grep -v ">" $initial_fasta )
do 
    access=$( timeout 60m python3 $additional_folder/access_py.py -s $sequence 2>>errors.log ) 
    exit_status=$?
    if [ -z "$access" ]  # If no value calculated, record NA
    then
        echo NA >> ./data/access.csv
    elif [[ "$access" == "THIS happened"* ]]   # If error occurred, record NA
    then
        echo NA >> ./data/access.csv
    elif [ "$exit_status" -eq "124" ]    # If calculation timed out, record NA
    then
        echo NA >> ./data/access.csv
    else
        echo $access >> ./data/access.csv
    fi
done

sed -i 's/nan/NA/g' ./data/access.csv  # make NA readable by R

rm -rf ./data/MFE.csv

############################################################################################################################

# Coding Potential - Ficket Score

############################################################################################################################

## (python script so cannot be added to $PATH)

$cpc2_directory/CPC2.py -i $initial_fasta -o ./data/output.txt >/dev/null 2>>errors.log    # Original version
#python3 $cpc2_directory/CPC2.py -i $initial_fasta -o output >/dev/null 2>>errors.log   # Uses updated biopython packages
cat ./data/output.txt | cut -f 4 > ./data/CPC2.csv
rm -rf ./data/output.txt


############################################################################################################################

# Obtain MSA siles RNAalifold scores, MFE, RNAcode score

###########################################################################################################################

var="$first_rna_id"
#rm -rf *.stk

while IFS=$'\t' read -r line
do
    ######## Obtain MSA files from multiz100way, unless they've already been downloaded
    rna_id='RNA'$var
    if [ -f ./data/maf/$rna_id.maf ]
    then
        cp ./data/maf/$rna_id.maf mafOut
    else            
        #until
        overBed="$line"  #CHECK IF WORKS WITH A VARIABLE AS INPUT INSTEAD OF A FILE!!   
	    echo "$mafFetch_exe hg38 multiz100way overBed mafOut 2>>errors.log ;" >> errors.log
	    $mafFetch_exe hg38 multiz100way overBed mafOut 2>>errors.log ;
	    #do sleep 4 ; done
        cp mafOut ./data/maf/$rna_id.maf 
    fi
	
    ######## If MSA contains data 
    file_len=$( cat mafOut | wc -l )
    if [ "$file_len" != "1" ]
    then
        
        ######## Run RNAcode
	    echo "$RNAcode_exe mafOut -o rnacode_output 2>>errors.log &> /dev/null" >>errors.log
	    $RNAcode_exe mafOut -o rnacode_output 2>>errors.log &> /dev/null
            
        # Takes the largest score, but records zero if no significant hits found.
        rnacode=$( grep -v "No significant coding regions found.\|#\|=" rnacode_output | grep . | tr -s ' ' | cut -d ' ' -f 10 | grep . | sort -V | tail -1 )
        [ -z $rnacode ] && rnacode=NA
	    

        ####### Convert MSA.maf files to STK format for RNAalifold and R-scape
        ./bin/maf_to_stockholm.sh $rna_id  # is it better to declare function have it on a different script??  

        ######## Generate secondary structure consensus sequence and associated MFE value
	    echo "$RNAalifold_exe -f S --noPS --aln RNA.stk >$rna_id.rnaalifold 2>>errors.log" >>errors.log
        #timeout 60m
	    $RNAalifold_exe -f S --noPS --aln RNA.stk > $rna_id'.rnaalifold' 2>>errors.log
        rna_score=$( cat $rna_id'.rnaalifold' | tail -1 | cut -d ' ' -f 2- | tr -d "(" | cut -d "=" -f 1 )

	    # Remove gap-only columns:
	    echo "esl-alimask -g --gapthresh 1 $rna_id.stk > maf/$rna_id.stk" >>errors.log
	    esl-alimask -g --gapthresh 1 $rna_id.stk > maf/$rna_id.stk
	    
        [ -z "$rna_score" ] && rna_score=NA
        echo $rnacode,$rna_score | tr -d ' ' >> rnacode.csv
	    
        exit_status=$?
	    
        ######## Run R-scape only if previous analyses didn't time out
        if [ "$exit_status" -ne "124" ]
        then
		    echo "$rscape_exe -E 100 -s $rna_id.stk >/dev/null 2>>errors.log" >>errors.log
            $rscape_exe -E 100 -s $rna_id.stk >/dev/null 2>>errors.log
        else
            touch $rna_id.cov
        fi
    
    else
        ######## Generate blanks if no MAF available
        touch $rna_id.cov
        echo NA,NA >> rnacode.csv
    fi

    (( var++ ))
    mv *.cov ./data/rscapedata &> /dev/null

done < ./data/coordinates


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
            file=./rscapedata/RNA$count_file.cov
            max=$(  grep -r 'GTp' $file       | cut -d '[' -f 2 | tr ']' ',' | cut -d ',' -f 2)  # Max covariance observed
            min=$(  grep -v "#" $file      | cut -f 4 | sort -n | head -1 | tr -d '-' )  # Min covariance observed
            eval=$( grep -m 1 "$min" $file | cut -f 5 )  # E-val of min covariance
            [[ "$eval" == "no significant pairs" ]] && eval=NA   
            echo $max,$eval >> rscape-dataset.csv
        else
            ######## If no covariance calculated/no data available
            echo NA,NA >> rscape-dataset.csv
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