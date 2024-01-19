# Function to convert MAF to Stockholm

maf_to_stockholm() {

######## Extract all the multiz100way alignment information 
grep -v "#\|score" data/maf/mafOut > data/maf/mafOut.txt

######## Split each file into individual alignments
awk -v RS= '{print > ("data/maf/mafOut-" NR ".txt")}' data/maf/mafOut.txt


max=$( ls data/maf/mafOut-*.txt | wc -l )
count=1

while [ $count -le $max ]
do
    ######## Parse through each file, in order 
    while IFS=$'\n' read -r line; do

        start=$( echo $line | tr -s ' ' | cut -d ' ' -f 1 )
        ######## If the line contains sequence information
        if (( $(echo "$start == s" | bc -l ) ))
        then
            ######## If this is the first file, then also extract sequence name
            if (( $(echo "$count == 1" | bc -l) )); then
                id=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | tr '.' '/' )
                seq=$(  echo $line | tr -s ' ' | cut -d ' ' -f 7 )
                
                if [ ! -z "$seq" ]; then
                    ######## Ensures all the sequences start in the same position in the STK file
                    printf "%-40s\n" $id >> data/maf/mafOut-0-seq
                    echo $seq >> data/maf/mafOut-$count-seq
                fi    

            else
                ######## Only extract sequence information in subsequent files
                seq=$( echo $line | tr -s ' ' | cut -d ' ' -f 7 )
                if [ ! -z "$seq" ]; then
                    echo $seq >> data/maf/mafOut-$count-seq
                fi
            fi
        else
            :
        fi

    done < data/maf/mafOut-$count.txt

    (( count++ ))  
done


######## Format STK file
if [ ! -d data/STK ]; then mkdir data/STK; fi
echo '# STOCKHOLM 1.0' > data/STK/$rna_id.stk
echo "#=GF ID $rna_id" >> data/STK/$rna_id.stk

# Concatenate sequences into a single file
files=()
while IFS= read -r file; do
    files+=("$file")
done < <(ls data/maf/mafOut-*-seq | sort -V)

paste -d'\0' "${files[@]}" > data/alignment

# Find the length of the sequence of the full sequence to remove partial alignments
human_len=$( grep "hg38" data/alignment | tr -s ' ' | cut -d ' ' -f 2 | wc -m )

while IFS=$'\n' read -r line; do     

    seq_len=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | wc -m )
    
    if [[ $seq_len -eq $human_len ]]; then                          # Only include full sequences          
        echo $line >> data/STK/$rna_id.stk
    fi

done < data/alignment


rm -rf data/maf/mafOut-*.txt
rm -rf data/maf/mafOut-*-seq
rm -rf data/alignment

}