negative_control() {

    local negative-control-position=$1

    chromo=$( echo $chr | tr -d "chr" )
    #Generate null sequences upstream of sequence
    left_end=$(( $start - "$negative-control-position" ))
    left_start=$(( $left_end - $left_length ))                   # length is based on exons to match selected functional exons
    
    if [[ $left_end -lt 0 ]] || [[ $left_start -lt 0 ]] ; then left_sequence= ; else left_sequence=$( grep -w "chromosome $chromo" $genome_csv | cut -f 2 | cut -c$left_start-$left_end ) ; fi
   
    if [ -z $left_sequence ]   # If no sequence extracted, then remove.
    then
        :
    elif [[ $left_sequence == *"N"* ]]   # If sequence extracted contains unknown nucleotides (N), then remove.
    then
        :
    else
        echo $chr,$left_start,$left_end,$left_sequence >> ./data/left-negative-control.csv
    fi

    #Generate null sequences downstream of sequence
    right_start=$(( $end + "$negative-control-position" ))
    right_end=$(( $right_start + $right_length ))
    if [[ $right_start -lt 0 ]] || [[ $right_end -lt 0 ]] ; then right_sequence= ; else right_sequence=$( grep -w "chromosome $chromo" $genome_csv | cut -f 2 | cut -c$right_start-$right_end ) ; fi
    
    if [ -z $right_sequence ]     # If no sequence extracted, then remove.
    then
        :
    elif [[ $right_sequence == *"N"* ]]   # If sequence extracted contains unknown nucleotides (N), then remove.
    then
        :
    else
        echo $chr,$right_start,$right_end,$right_sequence >> ./data/right-negative-control.csv
    fi

}

ncrna_bed= 
uniprot_bed=
###########################################################################################################################

# Create negative control protein-coding sequences from human genome

###########################################################################################################################

######## Generate sequences upstream and downstream of functional gene

negative_control_positions=(20000, 40000, 60000, 80000, 100000)

while IFS=, read chr start end left_len right_len   #start exon 1, end last exon, length exon 2 and 3
do
    for position in $[negative_control_positions[@]]; do
        negative_control "$position"
    done

done < coords-for-negative-control.csv

###########################################################################################################################

# Filter upstream negative control sequences using UniProt and GENCODE

###########################################################################################################################

######## Reformat data for bedtools
grep -v "Start" ./data/left-negative-control.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > ./data/left-coordinates.bed

######## Uniprot filtering
[ -f $dep_folder/bedtools ] && $dep_folder/bedtools intersect -a $uniprot_bed -b ./data/left-coordinates.bed > ./data/left-intersect.bed || bedtools intersect -a $uniprot_bed -b ./data/left-coordinates.bed > ./data/left-intersect.bed

######## Extract chromosome coordinates of overlapping negative control sequences
cat ./data/left-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > ./data/left-intersect
left_count=$( cat ./data/left-intersect | wc -l )

######## GENCODE filtering
[ -f $dep_folder/bedtools ] && $dep_folder/bedtools intersect -a $ncrna_bed -b ./data/left-coordinates.bed > ./data/left-intersect.bed || bedtools intersect -a $ncrna_bed -b ./data/left-coordinates.bed > ./data/left-intersect.bed

######## Extract chromosome coordinates of overlapping negative control sequences
cat ./data/left-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > ./data/left-intersect

######## Remove negative control sequences that overlap with known functional sequences
echo ID,Functional,Chromosome,Start,End,Sequence > ./data/negative-control-dataset.csv
count=$( tail -1 ./data/protein-exon2-dataset.csv | cut -d ',' -f 1 | tr -d "RNA" )

for line in $( cat ./data/left-negative-control.csv )
do
    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates ./data/left-intersect               # If overlap with protein-coding gene
    then
        :
    elif grep -q $coordinates ./data/left-intersect            # If overlap with ncRNA gene
    then
        :
    else
        (( count++ ))
        echo RNA$count,No,$line >> ./data/negative-control-dataset.csv
    fi
done

###########################################################################################################################

# Filter downstream negative control sequences using UniProt and GENCODE

###########################################################################################################################

######## Reformat data for bedtools
grep -v "Start" ./data/right-negative-control.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > ./data/right-coordinates.bed

######## Uniprot filtering
[ -f $dep_folder/bedtools ] && $dep_folder/bedtools intersect -a $uniprot_bed -b ./data/right-coordinates.bed > ./data/right-intersect.bed || bedtools intersect -a $uniprot_bed -b $d-right-coordinates.bed > $d-right-intersect.bed

######## Extract chromosome coordinates of overlapping negative control sequences
cat ./data/right-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > ./data/right-intersect
right_count=$( cat ./data/right-intersect | wc -l )

######## GENCODE filtering
[ -f $dep_folder/bedtools ] && $dep_folder/bedtools intersect -a $ncrna_bed -b ./data/right-coordinates.bed > ./data/right-intersect.bed || bedtools intersect -a $ncrna_bed -b ./data/right-coordinates.bed > ./data/right-intersect.bed

######## Extract chromosome coordinates of overlapping negative control sequences
cat ./data/right-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > ./data/right-intersect

for line in $( cat ./data/right-negative-control.csv )
do
    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates ./data/right-intersect      # If overlap with protein-coding gene
    then
        :
    elif grep -q $coordinates ./data/right-intersect       # If overlap with ncRNA gene
    then
        :
    else
        (( count++ ))
        echo RNA$count,No,$line >> ./data/negative-control-dataset.csv
    fi
done

######## Generate negative control FASTA file
grep -v "Start" ./data/negative-control-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > ./data/negative-control-seq.fa

###########################################################################################################################