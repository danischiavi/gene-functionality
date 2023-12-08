# Function to convert MAF to Stockholm

maf_to_stockholm() {

######## Extract all the multiz100way alignment information 
grep -v "#\|score" mafOut > mafOut.txt

######## Split each file into individual alignments
awk -v RS= '{print > ("mafOut-" NR ".txt")}' mafOut.txt

IFS=$'\n'
max=$( ls mafOut-*.txt | wc -l )
count=1

while [ $count -le $max ]
do
    ######## Parse through each file, in order 
    for line in $( cat mafOut-$count.txt )
    do
        start=$( echo $line | tr -s ' ' | cut -d ' ' -f 1 )
        ######## If the line contains sequence information
        if (( $(echo "$start == s" | bc -l ) ))
        then
            ######## If this is the first file, then also extract sequence name
            if (( $(echo "$count == 1" | bc -l) ))
            then
                name=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | tr '.' '/' )
                seq=$( echo $line | tr -s ' ' | cut -d ' ' -f 7 )
                if [ -z "$seq" ]
                then
                    :
                else
                    ######## Ensures all the sequences start in the same position in the STK file
                    printf "%-40s\n" $name >> mafOut-0-seq
                    echo $seq >> mafOut-$count-seq
                fi      
            else
                ######## Only extract sequence information in subsequent files
                seq=$( echo $line | tr -s ' ' | cut -d ' ' -f 7 )
                if [ -z "$seq" ]
                then
                    :
                else
                    echo $seq >> mafOut-$count-seq
                fi
            fi
        else
            :
        fi

    done
    count=$(( $count + 1 ))  

done

######## Format STK file
echo "# STOCKHOLM 1.0" > RNA.stk
echo >> RNA.stk
echo "#=GF ID $1" >> RNA.stk
echo >> RNA.stk

list=$( ls mafOut-*-seq | sort -V )
paste -d'\0' $list >> alignment

######## Find the length of the sequence of the full sequence
human_len=$( grep "hg38" alignment | tr -s ' ' | cut -d ' ' -f 2 | wc -m )

for line in $( cat alignment )
do

    seq_len=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | wc -m )
    if [[ $seq_len -eq $human_len ]]     # Only include full sequences and removes partial alignments
    then
        echo $line >> RNA.stk
    else
        :
    fi
done


rm -rf mafOut-*.txt
rm -rf mafOut-*-seq
rm -rf alignment

}