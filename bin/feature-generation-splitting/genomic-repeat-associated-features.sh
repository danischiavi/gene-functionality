#!/bin/bash
#
# Script Name: genomic-repeat-associated-features.sh
#
#
# Description: This script calculates: 
#                       - minimum and average distance to repeat 
#                       - closest non overlapping repetitive element 
#                       - genomic copy number  
# for protein-coding exons, lncRNA exons, short ncRNAs and negative control sequences. 
#
# Input: $1 is the dataset 
#        $2 is the fasta file for the sequences
#        $3 is the file which contains the sequence coordinates sorted and bed formatted
#        $4 is the location of the folder with the local databases/datasets or version specific executables
#       
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.

############################################################################################################################

initial_data=$1
initial_fasta=$2
sorted-coordinates.bed=$3
additional_folder=$4

####Checks 

## Local database for blastn

blast_total=$( ls $additional_folder/human_genome* | wc -l )
if [[ $blast_total -ge 6 ]]
then
    human_genome=$additional_folder/human_genome
else
    until [ -f $human_genome.nhr ] ; do read -p "Please enter custom name of human genome blastn database: " h ; human_genome=$additional_folder/$h ; echo ; done
fi

## blastn

if find $additional_folder/ -executable -type f | grep -q blastn
then
    blastn_exe=$additional_folder/blastn
elif command -v blastn &> /dev/null
then
    blastn_exe=blastn
else
    echo Please check that blastn is installed.
    exit 1
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


######## Local bedfile of Dfam hits for bedtools

if [ -f $additional_folder/dfam-hg38-sorted.bed ]        
then
    dfam_hits=$additional_folder/dfam-hg38-sorted.bed
else
    until [ -f $dfam_hits ] ; do read -p "Please enter custom name of Dfam non-redundant hits (bed): " d ; dfam_hits=$additional_folder/$d ; echo ; done
fi


#### File generation
echo Genome_copy_number,(Genome_complete_match)?, mmseqs_copy_number, nhmmer_copy_number > data/copy-number-$name.csv
echo Chromosome,Start,End,Dfam_min,Dfam_sum > data/dfam-distance-$name.csv

############################################################################################################################

# Calculate genomic copy number (identify number of matches within human genome)

############################################################################################################################

## Run blastn and MMseqs2

$blastn_exe -query $initial_fasta -db $human_genome -evalue 0.01 -out data/blastn-output.csv -outfmt "10 qaccver saccver pident" >/dev/null 2>>errors.log

mmseqs easy-search $initial_fasta data/raw/mmseqs/human_genome data/mmseqs-output data/tmp --search-type 3 --format-output query,target,evalue

## Run nhmmer and process Blastn and mmseqs results in order 

var=$first_rna_id                                                   # Global variables for given id number to first and last RNA sequence of the dataset
last_seq=$last_rna_id               

while [ $var -le $last_seq ]
do
    echo "$( grep -w -A 1 ">RNA$var" $initial_fasta )" > data/nhmmer-input.fasta
    nhmmer --tblout data/nhmmer-output -E 0.01 --noali data/nhmmer-input.fasta data/raw/GRCh38_p14_genomic.fna > /dev/null 2>&1 # -E <x> : report sequences <= this E-value threshold in output; don't output alignments
    nhmmer_hits=$( grep -v "#" data/nhmmer-output | wc -l )

    #minimap2 data/raw/GRCh38_p14_genomic.fna data/nhmmer-input.fasta >> data/minimap-output
    
    ######## Process each RNA in the file (by ID) as the counter increases
    total_blastn=$( grep -w "RNA$var" data/blastn-output.csv | wc -l)  
    total_mmseqs=$( grep -w "RNA$var" data/mmseqs-output | wc -l)  

    if [ -z "$total_blastn" ]; then total_blastn='NA'; else :; fi 
    if [ -z "$total_mmseqs" ]; then total_mmseqs='NA'; else :; fi 
    if [ -z "$nhmmer_hits" ]; then nhmmer_hits='NA'; else :; fi
    
    echo "$total_blastn, $total_mmseqs, $nhmmer_hits" >> data/copy-number-$name.csv
    (( var++ ))
   
done

## Telomere-to-telomere assembly

var=$first_rna_id                                                   # Global variables for given id number to first and last RNA sequence of the dataset
last_seq=$last_rna_id 
T2T_genome='/Volumes/archive/userdata/student_users/danielaschiavinato/dani-scratch/features-of-functional-human-genes/data/raw/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna'              

while [ $var -le $last_seq ]
do
    echo "$( grep -w -A 1 ">RNA$var" $initial_fasta )" > data/nhmmer-input.fasta
    nhmmer --tblout data/nhmmer-output -E 0.01 --noali data/nhmmer-input.fasta $T2T_genome > /dev/null 2>&1 # -E <x> : report sequences <= this E-value threshold in output; don't output alignments
    nhmmer_hits_T2T=$( grep -v "#" data/nhmmer-output | wc -l )
    
    if [ -z "$nhmmer_hits_T2T" ]; then nhmmer_hits_T2T='NA'; else :; fi

    echo "$nhmmer_hits_T2T" >> data/T2T-copy-number-$name.csv

    (( var++ ))

done



############################################################################################################################

# Distance to nearest transposable element (Dfam)

############################################################################################################################

######## Format data for bedtools
grep -v 'Chromosome' $initial_data |  cut -f 3,4,5 -d "," | tr ',' '\t' > data/initial-data-$name.bed
sort -k1,1 -k2,2n data/initial-data.bed | tr ' ' '\t' > data/initial-data-sorted-$name.bed 

######## Upstrean hits
$bedtools_exe closest -a data/initial-data-sorted-$name.bed -b $dfam_hits -io -D ref -iu > data/dfam-downstream-$name.bed 2>>errors.log

######## Downstream hits
$bedtools_exe closest -a data/initial-data-sorted-$name.bed -b $dfam_hits -io -D ref -id > data/dfam-upstream-$name.bed 2>>errors.log

paste <( cut -f 1,2,3,7 data/dfam-downstream-$name.bed ) <( cut -f 7 data/dfam-upstream-$name.bed ) --delimiters '\t' > data/dfam-combined-$name.bed


######## Calculate sum of upstream and downstream, and combine into one file
while IFS=$'\t' read -r chr start end downstream upstream
do
    chr_numb=$( echo "$chr" | cut -c 4- )
    upstream=$( echo "$upstream" | tr -d '-' )
    [ -z "$downstream" ] && downstream=0                                                # Implies that no region up or downstream, rather than data missing.
    [ -z "$upstream" ] && upstream=0 
    sum=$(( $downstream+$upstream ))
    if [[ "$upstream" -lt "$downstream" ]]; then                                              # Taking into account forward and reverse strand
        echo $chr,$start,$end,$upstream,$sum >> data/dfam-distance-$name.csv
    else
        echo $chr,$start,$end,$downstream,$sum >> data/dfam-distance-$name.csv    
    fi

done < data/dfam-combined-$name.bed

#### Remove unnessesary files
mv dfam_downstream.bed additional-output/
mv dfam_upstream.bed additional-output/
rm -rf data/dfam-combined-$name.bed
rm -rf data/dfam-upstream-$name.bed
rm -rf data/dfam-downstream-$name.bed
rm -rf data/initial-data-$name.bed

rm -rf data/nhmmer-input.fasta
rm -rf data/nhmmer-output





######## Function for matching up dfam-distance to ncRNA in R (bedtools is unordered)

reformat_dfam() {

R --save << RSCRIPT
df1 <- read.csv("file1.csv", stringsAsFactors=F)
df2 <- read.csv("file2.csv", stringsAsFactors=F)
for(i in 1:nrow(df2)){
        index <- grep(df2[i, 2], df1[,4])
        df1[index,'snp_num'] <- df2[i,4]
        df1[index,'snp_ave'] <- df2[i,5]
}
write.csv(df1, "file3.csv", quote=F, row.names=F)
RSCRIPT

}

######## Reorder dfam-distance data
cat $initial_data > file1.csv
cat dfam-distance.csv > file2.csv
reformat_dfam >/dev/null 2>>errors.log
cat file3.csv | cut -d ',' -f 7,8 > ncrna-dfam-distance.csv

rm -rf dfam-distance.csv
rm -rf rnacentraloutput_FINAL_sorted1.bed
rm -rf rnacentraloutput_FINAL.bed
rm -rf rnacentraloutput_FINAL_sorted.bed
rm -rf file1.csv
rm -rf file2.csv
rm -rf file3.csv

mv snp-intersection.bed additional-output/
mv snp-intersection-average.csv additional-output/