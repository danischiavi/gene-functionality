# Function to convert MAF to Stockholm

maf_to_stockholm() {

    local input_file=$1
    local file_id=$2

    mkdir -p data/maf-to-stk
    local temp_dir='data/maf-to-stk'
   
    ######## Extract all the multiz100way alignment information 
    grep -v "#\|score" $input_file > $temp_dir/aligns.txt

    ######## Split each file into individual alignments
    awk -v RS= '{print > ("'"$temp_dir"'/align-" NR ".txt")}' "$temp_dir"/aligns.txt


    max=$( ls $temp_dir/align-*.txt | wc -l )
    count=1

    while [ $count -le $max ]; do
        ######## Parse through each file, in order 
        while IFS=$'\n' read -r line; do

            start=$( echo "$line" | tr -s ' ' | cut -d ' ' -f 1 )
            ######## If the line contains sequence information
            if (( $(echo "$start == s" | bc -l ) )); then
                ######## If this is the first file, then also extract sequence name
                if (( $(echo "$count == 1" | bc -l) )); then
                
                    id=$(awk '{print $2"/"$3"-"($3 + $4)}' <<< "$line")
                    seq=$(awk '{print $7}' <<< "$line")
                
                    if [ ! -z "$seq" ]; then
                        ######## Ensures all the sequences start in the same position in the STK file
                        printf "%-90s\n" $id >> $temp_dir/align-0-seq
                        echo $seq >> $temp_dir/align-$count-seq
                
                    fi    

                else
                ######## Only extract sequence information in subsequent files
                    seq=$(awk '{print $7}' <<< "$line")
                
                    if [ ! -z "$seq" ]; then echo $seq >> $temp_dir/align-$count-seq; fi
        
                fi
            fi

        done < $temp_dir/align-$count.txt

        (( count++ ))  

    done


    ######## Format STK file

    if [ ! -d data/STK/$name ]; then mkdir data/STK/$name; fi

    local stk_file=data/STK/$name/$file_id.stk
    
    echo '# STOCKHOLM 1.0' > $stk_file
    echo "#=GF ID $file_id" >> $stk_file

    # Concatenate sequences into a single file
    files=()

    while IFS=$'\n' read -r file; do

        files+=("$file")

    done < <(ls $temp_dir/align-*-seq | sort -V)


    paste -d'\0' "${files[@]}" > $temp_dir/alignment

    # Find the length of the sequence of the full sequence to remove partial alignments
    human_len=$( grep "hg38" $temp_dir/alignment| tr -s ' ' | cut -d ' ' -f 2 | wc -m )

    while IFS=$'\n' read -r line; do     

        seq_len=$( echo "$line" | tr -s ' ' | cut -d ' ' -f 2 | wc -m )
    
        if [[ $seq_len -eq $human_len ]]; then echo $line >> $stk_file; fi  # Only include full sequences

    done < $temp_dir/alignment

    echo "//" >> $stk_file

    rm -rf $temp_dir

}