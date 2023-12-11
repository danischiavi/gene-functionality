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

if find $additional_folder/ -executable -type f | grep -q bigWigSummary
then
    bigWigSummary_exe=$additional_folder/bigWigSummary
elif command -v bigWigSummary &> /dev/null
then
    bigWigSummary_exe=bigWigSummary
else
    echo Please check that bigWigSummary is installed.
    exit 1
fi

## Local bigwig file of UCSC phyloP scores for BigWigSummary

if [ -f $additional_folder/hg38.phyloP100way.bw ] 
then
    phylo_bw=$additional_folder/hg38.phyloP100way.bw
else
    until [ -f $phylo_bw ] ; do read -p "Please enter custom name of hg38 conservation scoring by phyloP (bigWig): " phylo ; phylo_bw=$additional_folder/$phylo ; echo ; done
fi

##Local bigwig file of UCSC phastCons scores for BigWigSummary

if [ -f $additional_folder/hg38.phastCons100way.bw ] 
then
    phast_bw=$additional_folder/hg38.phastCons100way.bw 
else
    until [ -f $phast_bw ] ; do read -p "Please enter custom name of hg38 conservation scoring by phastCons (bigWig): " phast ; phast_bw=$additional_folder/$phast ; echo ; done
fi

##Local bigwig file of UCSC phastCons scores for BigWigSummary

if [ -f $additional_folder/hg38.phastCons100way.bw ] 
then
    zoonomia_phyloP_bw=./data/241-mammalian-2020v2.bigWig
else
    until [ -f $zoonomia_phyloP_bw ] ; do read -p "Please enter custom name of hg38 conservation scoring by phastCons (bigWig): " phast ; phast_bw=$additional_folder/$phast ; echo ; done
fi

##Local bigwig file of 111 mammal GERP scores for BigWigSummary

if [ -d $additional_folder/gerp-mammals-index/ ] && [ -f $additional_folder/gerp-mammals-index/gerp1.bed.gz ] 
then
    mammals_bed=$additional_folder/gerp-mammals-index
else
    until [ -f $mammals_bed ] ; do read -p "Please enter custom name of sorted 111 mammal GERP scores (bedgraph): " mammal ; mammals_bed=
$additional_folder/$mammal ; echo ; done
fi

######## bedtools

if find $additional_folder/ -executable -type f | grep -q bedtools
then
    bedtools_exe=$additional_folder/bedtools
elif command -v bedtools &> /dev/null
then
    bedtools_exe=bedtools
else
    echo Please check that bedtools is installed.
    exit 1
fi

############################################################################################################################

# Calculating sequence conservation features

############################################################################################################################

echo MeanPhyloP-Z,MaxPhyloP-Z,MeanPhyloP,MaxPhyloP,MeanPhastCons,MaxPhastCons > ./data/$name_set-conservation.csv
echo mammals_mean_gerp,mammals_max_gerp > ./data/$name_set-gerp.csv

var="$first_rna_id"

while IFS=$'\t' read -r chr start end
do
    #echo '#####################RNA'$var'#####################' >> errors.log #D:i think I dont need the var in this section. Maybe another to append to error.log?
    
    
    ######## Obtain phastCons (pc) and phyloP (pp) values for each set of chromosome coordinates
	echo "mean_pp = $bigWigSummary_exe -type=mean $phylo_bw $chr $start $end 1 2>&1" >> errors.log
        mean_pp=$( $bigWigSummary_exe -type=mean $phylo_bw $chr $start $end 1 2>&1 )
        max_pp=$(   $bigWigSummary_exe -type=max $phylo_bw $chr $start $end 1 2>&1 )
        mean_pc=$( $bigWigSummary_exe -type=mean $phast_bw $chr $start $end 1 2>&1 )
        max_pc=$(   $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1 )
	echo "max_pc = $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1" >> errors.log
	
    ######## Extract PhyloP (pp) values from 241-way mammalian alignment for each set of chromosome coordinates
    zoonomia_mean_pp=$( $bigWigSummary_exe -type=mean $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
    zoonomia_max_pp=$(   $bigWigSummary_exe -type=max $zoonomia_phylo_bw $chr $start $end 1 2>&1 )
    

    ######## Convert any missing data to NAs (mean = mn, max = mx, phyloP = p, phastCons = c)
        test_mnp=$( echo $mean_pp | wc -w )
        if [[ "$test_mnp" -eq "1" ]] ; then : ; else mean_pp=NA ; fi
	
        test_mxp=$( echo $max_pp | wc -w )
        if [[ "$test_mxp" -eq "1" ]] ; then : ; else max_pp=NA ; fi
	
        test_mnc=$( echo $mean_pc | wc -w )
        if [[ "$test_mnc" -eq "1" ]] ; then : ; else mean_pc=NA ; fi
	
        test_mxc=$( echo $max_pc | wc -w )
        if [[ "$test_mxc" -eq "1" ]] ; then : ; else max_pc=NA ; fi
	
        test_mnc=$( echo $zoonomia_mean_pp | wc -w )
        if [[ "$test_mnc" -eq "1" ]] ; then : ; else mean_pc=NA ; fi
	
        test_mxc=$( echo $zoonomia_max_pp | wc -w )
        if [[ "$test_mxc" -eq "1" ]] ; then : ; else max_pc=NA ; fi

    echo $zoonomia_mean_pp,$zoonomia_max_pp$mean_pp,$max_pp,$mean_pc,$max_pc >> ./data/$name_dataset-conservation.csv

    
    
    
    ######## Obtain GERP Scores for each set of chromosome coordinates
    ## Format data for bedtools
    #need to strip "chr" from input.bed & from "gerp$chr.bed.gz"
    chromo=$( echo $chr | tr -d 'chr' )               
    echo -e "$chromo\t$start\t$end" > ./data/input-gerp.bed   # Input for bedtools (GERP) -e option for tab recognition
	echo "[chr:$chr] [chromo:$chromo] [start:$start] [end:$end]" >> errors.log 

	echo "$bedtools_exe map -a ./data/input-gerp.bed -b $mammals_bed/gerp$chromo.bed.gz -c 4,4 -o mean,max -sorted 2>>/dev/null" >> errors.log  #  map: Apply a function to a column for each overlapping interval
	mammal_output=$( $bedtools_exe map -a ./data/input-gerp.bed -b $mammals_bed/gerp$chromo.bed.gz -c 4,4 -o mean,max -sorted 2>> errors.log)     # -c Specify the column from the B file to map onto intervals in A 
    mean_mammal=$( echo $mammal_output | cut -f 4 )   # Mean GERP Score
    max_mammal=$( echo $mammal_output | cut -f 5 )    # Max GERP Score
	
    ######## Convert any missing data to NAs
    if [[ ! "$mean_mammal" == "." ]] ; then : ; else mean_mammal=NA ; fi
    if [[ ! "$max_mammal" == "." ]] ; then : ; else max_mammal=NA ; fi
	
    echo $mean_mammal,$max_mammal >> ./data/gerp.csv
 
[[[]	
	mkdir -p aliout
        var=$(($var+1))
        rm -rf *.sorted.cov
        mv *.cov rscapedata &> /dev/null
        rm -rf overBed
        rm -rf blank.txt
        #rm -rf *.pdf
        #rm -rf *ss.ps
        #rm -rf *.svg
        #rm -rf *.surv
        #rm -rf *.sto
	mv RNA*.ps  aliout/
	mv RNA*.stk aliout/
	mv RNA*svg  aliout/
	mv RNA*sto  aliout/
	mv RNA*pdf  aliout/
	mv RNA*surv  aliout/
	mv *rnaalifold aliout/
        rm -rf *.power
        rm -rf rnacode_output

done < ./data/coordinates

zip -r maf maf &> /dev/null
rm -rf maf/