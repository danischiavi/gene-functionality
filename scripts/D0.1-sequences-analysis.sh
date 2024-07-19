gencode_annotations=/Volumes/archive/userdata/student_users/danielaschiavinato/dani-scratch/gene-functionality/data/raw/gencode-lncRNA-annotations.gff3

dir=data/gencode-lncRNA
mkdir -p "$dir"

# Gene list 
awk -F'\t' '$3 == "gene" {split($NF, arr, ";"); for (i in arr) {if (match(arr[i], /^ID=([^;]+)/, id)) {print id[1]}}}' "$gencode_annotations" > data/gencode-lncRNA/gencode-lncrna-list 	# maybe add it to readmedata 

#### Temporary files ####
exon_one_info="${dir}/gencode-lncrna-exon1-info"
exon_two_info="${dir}/gencode-lncrna-exon2-info"
exon_one_tmp="${dir}/gencode-lncrna-exon1-tmp"
exon_two_tmp="${dir}/gencode-lncrna-exon2-tmp"
combined_seq="${dir}/combined-seq"

gff2Info() {                                                                                

    local exons=$1
	local id=$2

    coords_one=$( awk 'NR==1 {print $1, $4, $5, $7}' "$exons" )                                                  
    coords_two=$( awk 'NR==2 {print $1, $4, $5}' "$exons" )                                                    
      												
    IFS=' ' read -r chr start_one end_one strand <<< "$coords_one"                                                 # Starting and end coordinates of exons
    IFS=' ' read -r _ start_two end_two <<< "$coords_two"
	start_one_bed=$(( start_one - 1 ))															# to follow bed format required for bedtools getfasta 
	start_two_bed=$(( start_two - 1 ))
                                                                       
    len_one=$(( end_one - start_one_bed ))
	len_two=$(( end_two - start_two_bed ))           
       
	# Exclude sequences out of length limits 
	if { [ "$len_two" -lt "$lower_limit" ] || [ "$len_one" -lt "$lower_limit" ];}; then

		small_size+=("$id")

    elif { [ "$len_two" -gt "$upper_limit" ] && [ "$len_one" -gt "$upper_limit" ]; }; then 

		large_size+=("$id")	
	
	else 

		echo -e "${chr}\t${start_one}\t${end_one}\t${id}\t.\t${strand}" >> "$exon_one_info"		# following bed format 
        echo -e "${chr}\t${start_two}\t${end_two}\t${id}\t.\t${strand}" >> "$exon_two_info" 

    fi 
}



ambiguity() {

    local input_file=$1

    while IFS=$'\t' read -r name sequence; do

		id="${name%%::*}"
        # Calculate the ambiguity percent
        total_nucleotides=${#sequence}
        amb_nucleotides=$(echo "$sequence" | grep -o -i '[RYMKSWBDHVN]' | wc -l)
        ambiguity_percent=$((amb_nucleotides * 100 / total_nucleotides))

        # Check if ambiguity is higher than 5
        if [ "$ambiguity_percent" -gt 5 ]; then
            high_ambiguity_ids+="$id"
		else
			leftOver+="$id "
        fi
    done < "$input_file"

    echo "ambiguous: ${#high_ambiguity_ids[@]}"
	echo "leftOver: ${#leftOver[@]}"
}

reformat_file() {

    local input_file=$1
    local output_file=$2

    awk -F'\t' '
    function calc_ambiguity(seq,   i, amb_nucleotides, total_nucleotides, ambiguity_percent) {
        amb_nucleotides = 0
        total_nucleotides = length(seq)
        
        for (i = 1; i <= total_nucleotides; i++) {
            if (index("RYMKSWBDHVNrymkswbdhvn", substr(seq, i, 1))) {
                amb_nucleotides++
            }
        }
        
        ambiguity_percent = (amb_nucleotides * 100) / total_nucleotides
        return ambiguity_percent
    	
		}

    {
        split($1, parts, "::")
        split(parts[2], coords, ":")
        gene_id = parts[1]
        chromosome = coords[1]
        split(coords[2], range, "-")
        start = range[1] + 1 
        end = substr(range[2], 1, index(range[2], "(") -1)

        sequence = $2
        ambiguity = calc_ambiguity(sequence)

        if (ambiguity <= 5) {
            printf("NA,%s,%s,%s,%s,%s\n", chromosome, start, end, sequence, gene_id)
        }
    }' "$input_file" >> "$output_file"

	# start + 1 since is following bed format 

    rm -rf "$input_file"
}


declare -a "IDs_gencode_lncrna=()"
mapfile -t IDs_gencode_lncrna < data/gencode-lncRNA/gencode-lncrna-list

  
# Get GeneIDs that has already been used to train/test predicting model 
declare -a "selected_IDs=()"
values=$(awk -F, '{print $7}' data/datasets/gencode-lncRNA/gencode-lncrna-exon1-dataset.csv)
IFS=$'\n' read -r -d '' -a selected_IDs <<< "$values"

declare -a "excluded_IDs=()"
for id in "${IDs_gencode_lncrna[@]}"; do

	if [[ ! " ${selected_IDs[@]} " =~ " $id " ]]; then  

		excluded_IDs+=("$id")
	
	fi

done

declare -a "number_exons=()" 
declare -a "small_size=()"
declare -a "large_size=()" 
declare -a "high_ambiguity_ids=()" 
declare -a "leftOver=()"

for id in "${excluded_IDs[@]}"; do

	awk -F'\t' -v id="$id" '$0 ~ "gene_id="id && $3 == "exon"' "$gencode_annotations" > data/gencode-lncRNA/gencode-exons

    if [ "$(wc -l < data/gencode-lncRNA/gencode-exons)" -ge 2 ]; then
	
		gff2Info data/gencode-lncRNA/gencode-exons "$id"

	else

		number_exons+=("$id")

	fi    
		
done

echo "${#number_exons[@]}" 
echo "${#small_size[@]}" 
echo "${#large_size[@]}"

# Extract sequences from genome ##
bedtools getfasta -fi "$genome_seq" -bed "$exon_one_info" -fo "$exon_one_tmp" -s -name -tab  
			# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
			# -name Use the name field and coordinates for the FASTA header
			# tab Report extract sequences in a tab-delimited format instead of in FASTA format.

bedtools getfasta -fi "$genome_seq" -bed "$exon_two_info" -fo "$exon_two_tmp" -s -name -tab 


awk 'NR==FNR{a[NR]=$2; next} {print $1 "\t" $2 a[FNR]}' "$exon_one_tmp" "$exon_two_tmp" > "$combined_seq"

ambiguity "$combined_seq"


######## PLOT IN PYTHON ##########
import matplotlib.pyplot as plt

labels = ['< 2 exons', 'exons > 1500nt', 'exons < 74nt', 'seq ambiguity']
sizes = [4191, 60, 3937, 0]
selected_ids = 11971
excluded_total = sum(sizes)
sizes = [20198 - excluded_total] + sizes  # Include remaining lncRNAs that are not excluded
labels = ['Selected IDs'] + labels

# Colors
colors = ['green', 'pink', 'red', 'orange', 'yellow']

# Function to display percentages and raw values
def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return f'{pct:.1f}%\n({val:,})'
    return my_format

# Plotting
fig, ax = plt.subplots(figsize=(7, 7))
wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors, autopct=autopct_format(sizes), 
                                  startangle=140, pctdistance=0.6, labeldistance=0.8)

# Adjusting the position of labels to avoid overlap
for i, text in enumerate(texts):
    text.set_horizontalalignment('center')
    text.set_verticalalignment('center')

# Adding a title
plt.title('Gencode lncRNA Genes Filtering Summary')

# Equal aspect ratio ensures that pie is drawn as a circle
plt.axis('equal')

# Save the plot
plt.savefig('lncRNA_genes_filtering_summary.png')

# Display the plot
plt.show()