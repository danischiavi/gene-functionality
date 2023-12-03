#!/bin/bash
#
# Script Name: protein-and-rna-specific-features.sh
#
#
# Description: This script runs all the scripts to calculates the chosen features for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences. 
#
# Input: $1 is the file identifier for the dataset and fasta file. Eg: 200702-functional-ncrna
#        $2 is the location of the folder with the local databases/datasets or version specific executables
# 
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.

############################################################################################################################

initial_fasta= 
additional_folder=

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

#### File declaration
echo InteractionMIN,InteractionAVE > interaction-intermediate.csv
echo MFE > MFE.csv
echo Accessibility > access.csv
echo RNAcode_score,RNAalifold_score > rnacode.csv

############################################################################################################################

# RNA:RNA interactions

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

grep "(" rnafold-output | rev | cut -d "(" -f 1 | rev | tr -d ")" | tr -d " " >> MFE.csv
awk '!NF{$0="0"}1' MFE.csv > MFE-final.csv
rm -rf *.ps
rm -fr MFE.csv 


## Accessibility calculation
for sequence in $( grep -v ">" $initial_fasta )
do 
    access=$( timeout 60m python3 $additional_folder/access_py.py -s $sequence 2>>errors.log ) 
    exit_status=$?
    if [ -z "$access" ]  # If no value calculated, record NA
    then
        echo NA >> access.csv
    elif [[ "$access" == "THIS happened"* ]]   # If error occurred, record NA
    then
        echo NA >> access.csv
    elif [ "$exit_status" -eq "124" ]    # If calculation timed out, record NA
    then
        echo NA >> access.csv
    else
        echo $access >> access.csv
    fi
done

sed -i 's/nan/NA/g' access.csv  # make NA readable by R


############################################################################################################################

# Coding Potential

############################################################################################################################

## Fickett score calculation (python script so cannot be added to $PATH)

$cpc2_directory/CPC2.py -i $initial_fasta -o output.txt >/dev/null 2>>errors.log    # Original version
#python3 $cpc2_directory/CPC2.py -i $initial_fasta -o output >/dev/null 2>>errors.log   # Uses updated biopython packages
cat output.txt | cut -f 4 > CPC2.csv
rm -rf output.txt





###TO BE CONTINUED###
######## Create folder for generated r-scape output
rm -rf rscapedata && mkdir -p rscapedata

######## If multiz100way files already downloaded, use them instead of re-downloading 
if [ -f maf.zip ] ; then unzip maf.zip &> /dev/null && rm maf.zip ; else mkdir -p maf ; fi

######## Reformats chromosome coordinates for mafFetch
rm -rf coordinates
grep -v 'Chromosome' $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' > coordinates

######## Obtain multiple alignment files and calculate covariance
var=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
#rm -rf *.stk





while IFS=$'\t' read -r line
do
    ######## Obtain MSA files from multiz100way, unless they've already been downloaded
        rna_id='RNA'$var
        if [ -f maf/$rna_id.maf ]
        then
            cp maf/$rna_id.maf mafOut
        else            
            #until
	    echo "$mafFetch_exe hg38 multiz100way overBed mafOut 2>>errors.log ;" >> errors.log
	    $mafFetch_exe hg38 multiz100way overBed mafOut 2>>errors.log ;
	    #do sleep 4 ; done
            cp mafOut maf/$rna_id.maf 
        fi
	
        file_len=$( cat mafOut | wc -l )
        
        ######## If MSA contains data, convert is to STK format
        if [ "$file_len" != "1" ]
        then
            maf_to_stockholm $rna_id 
	    
            ######## Run RNAcode using original mafOut file
	    echo "$RNAcode_exe mafOut -o rnacode_output 2>>errors.log &> /dev/null" >>errors.log
	    $RNAcode_exe mafOut -o rnacode_output 2>>errors.log &> /dev/null
            
            ######## Takes the largest score, but records zero if no significant hits found.
            rnacode=$( grep -v "No significant coding regions found.\|#\|=" rnacode_output | grep . | tr -s ' ' | cut -d ' ' -f 10 | grep . | sort -V | tail -1 )
            [ -z $rnacode ] && rnacode=NA
	    
            ######## Generate secondary structure consensus sequence and associated MFE value
	    echo "$RNAalifold_exe -q -f S --noPS --aln-stk=$rna_id RNA.stk >$rna_id.rnaalifold 2>>errors.log" >>errors.log
            #timeout 60m
	    $RNAalifold_exe -q -f S --noPS --aln-stk=$rna_id RNA.stk > $rna_id'.rnaalifold' 2>>errors.log
            rna_score=$( cat $rna_id'.rnaalifold' | tail -1 | cut -d ' ' -f 2- | tr -d "(" | cut -d "=" -f 1 )

	    #Remove gap-only columns:
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

done < coordinates