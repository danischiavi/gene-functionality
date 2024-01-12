#!/bin/bash
#
# Script Name: seq-conservation-features.sh
#
# Description: This script calculates the PhyloP, PhyloCons and GERP for for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences.
#
# Input: $1 is the dataset 
#        $2 is the location of the folder with the local databases/datasets or version specific executables
# 
# Any additional required files, directories or dependencies will be requested when the script is run and require manual
# input.
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.
#

###########################################################################################################################

####Checks 

## bigWigSummary

if find bin/ -executable -type f | grep -q bigWigSummary
then
    bigWigSummary_exe=bin/bigWigSummary
elif command -v bigWigSummary &> /dev/null
then
    bigWigSummary_exe=bigWigSummary
else
    echo Please check that bigWigSummary is installed.
    exit 1
fi

## Local bigwig file of UCSC 100way phyloP scores for BigWigSummary

if [ -f data/raw/hg38.phyloP100way.bw ] ## UPDATE
then
    phylo_bw=data/raw/hg38.phyloP100way.bw
else
    until [ -f $phylo_bw ] ; do read -p "Please enter custom name of hg38 conservation scoring by phyloP (bigWig): " phylo ; phylo_bw=$additional_folder/$phylo ; echo ; done
fi

##Local bigwig file of UCSC phastCons scores for BigWigSummary

if [ -f $additional_folder/hg38.phastCons100way.bw ] ### update
then
    phast_bw=$additional_folder/hg38.phastCons100way.bw 
else
    until [ -f $phast_bw ] ; do read -p "Please enter custom name of hg38 conservation scoring by phastCons (bigWig): " phast ; phast_bw=$additional_folder/$phast ; echo ; done
fi

##Local bigwig file of UCSC 241way phastP scores for BigWigSummary

if [ -f data/raw/241-mammalian-2020v2.bigWig ]  
then
    zoonomia_phyloP_bw=data/raw/241-mammalian-2020v2.bigWig
else
    until [ -f $zoonomia_phyloP_bw ] ; do read -p "Please enter custom name of hg38 conservation scoring by phastCons (bigWig): " phast ; phast_bw=$additional_folder/$phast ; echo ; done
fi

##Local bigwig file of 111 mammal GERP scores for BigWigSummary ## UPDATE

if [ -d data/raw/gerp-mammals-index/ ] && [ -f $additional_folder/gerp-mammals-index/gerp1.bed.gz ] 
then
    mammals_bed=$additional_folder/gerp-mammals-index
else
    until [ -f $mammals_bed ] ; do read -p "Please enter custom name of sorted 111 mammal GERP scores (bedgraph): " mammal ; mammals_bed=
$additional_folder/$mammal ; echo ; done
fi


############################################################################################################################

# Calculating sequence conservation features

############################################################################################################################

echo Zoonomia-MeanPhyloP,Zoonomia-MaxPhyloP,MeanPhyloP,MaxPhyloP,MeanPhastCons,MaxPhastCons > data/conservation-$name.csv
echo GERP_mean,GERP_max > data/gerp-$name.csv

while IFS=',' read -r _ _ chr start end _
do
    ######## Obtain phastCons (pc) and phyloP (pp) values for each set of chromosome coordinates
	echo "mean_pp = $bigWigSummary_exe -type=mean $phylo_bw $chr $start $end 1 2>&1" >> errors.log
        mean_pp=$( $bigWigSummary_exe -type=mean $phylo_bw $chr $start $end 1 2>&1 )
        max_pp=$(   $bigWigSummary_exe -type=max $phylo_bw $chr $start $end 1 2>&1 )
        mean_pc=$( $bigWigSummary_exe -type=mean $phast_bw $chr $start $end 1 2>&1 )
        max_pc=$(   $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1 )
	echo "max_pc = $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1" >> errors.log
	
    ######## Extract PhyloP (pp) values from 241-way mammalian alignment for each set of chromosome coordinates
    zoonomia_mean=$( $bigWigSummary_exe -type=mean $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
    zoonomia_max=$(   $bigWigSummary_exe -type=max $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
    
    ######## Obtain GERP Scores for each set of chromosome coordinates
    gerp_mean=$( $bigWigSummary_exe -type=mean $gerp_bw $chr $start $end 1 2>&1 )
    gerp_max=$(   $bigWigSummary_exe -type=max $gerp_bw $chr $start $end 1 2>&1 )

    ######## Convert any missing data to NAs (mean = mn, max = mx, phyloP = p, phastCons = c) --> D: declare a short function for this
        test_mnp=$( echo $mean_pp | wc -w )
            if [[ "$test_mnp" -eq "1" ]]; then : ; else mean_pp=NA ; fi
	
        test_mxp=$( echo $max_pp | wc -w )
        if [[ "$test_mxp" -eq "1" ]]; then : ; else max_pp=NA ; fi
	
        test_mnc=$( echo $mean_pc | wc -w )
        if [[ "$test_mnc" -eq "1" ]]; then : ; else mean_pc=NA ; fi
	
        test_mxc=$( echo $max_pc | wc -w )
        if [[ "$test_mxc" -eq "1" ]]; then : ; else max_pc=NA; fi
	
        test_mnc=$( echo $zoonomia_mean | wc -w )
        if [[ "$test_mnc" -eq "1" ]]; then : ; else mean_pc=NA ; fi
	
        test_mxc=$( echo $zoonomia_max | wc -w )
        if [[ "$test_mxc" -eq "1" ]]; then : ; else max_pc=NA ; fi

        test_mnc=$( echo $gerp_mean | wc -w )
        if [[ "$test_mnc" -eq "1" ]]; then : ; else gerp_mean=NA ; fi
	
        test_mxc=$( echo $gerp_max | wc -w )
        if [[ "$test_mxc" -eq "1" ]]; then : ; else gerp_max=NA ; fi
    
    echo $zoonomia_mean,$zoonomia_max,$mean_pp,$max_pp,$mean_pc,$max_pc >> data/conservation-$name.csv
    echo $gerp_mean,$gerp_max >> data/gerp-$name.csv
 	
done < "$initial_data"


