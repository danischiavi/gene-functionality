# OBTAIN RAW DATA FOR ANALYSIS: 

# HUMAN GENOME: 

    * Download and unzip GFF file for human genome. 15/12/2023
    cd data/raw
    wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz   
    gunzip GRCh38_latest_genomic.gff.gz

    * Download, unzip FNA and format file for human genome. 15/12/2023
    cd data/raw
    wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz   && 
    gunzip GRCh38_latest_genomic.fna.gz
    fasta_formatter -i data/raw/GRCh38_p14_genomic.fna -o data/raw/GRCh38_interim.csv -t 

    # Remove any scaffolds or alternative chromosomes
    grep "NC_" data/raw/GRCh38_interim.csv > data/raw/GRCh38.p14_genome.csv
    rm data/raw/GRCh38_interim.csv


# RETRIEVAL OF FUNCTIONAL GENES:
# PROTEIN-CODING-GENES    
* A text file for all protein-coding genes from HGNC (Braschi et al., 2019), called gene_with_protein_product.txt, was downloaded and processed as described below.

    # Download text file for all protein-coding genes from HGNC
    wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt

    # Pipes one and two removes sequences encoded in the Y and Mt chromosomes or are no longer approved genes.
    # Pipes three to five extracts the columns associated with RefSeq IDs, and removes sequences associated with more than one ID.

    grep -v "mitochondrially encoded\|Entry Withdrawn\|Yq" data/raw/gene_with_protein_product.txt | grep "Approved" | cut -f 24 | grep -v "|" | grep . > data/raw/hgnc-protein-coding-RefSeq.txt



# NON-CODING GENES:

* To obtain the chromosome coordinates for each RNAcentral ncRNA, the Homo sapiens GRCh38 bed file was obtained from the RNAcentral FTP directory, which contains the chromosome coordinates for all human ncRNAs (The RNAcentral Consortium, 2019). 19/12/2024

wget -P data/raw ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/22.0/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz
gunzip -c data/raw/homo_sapiens.GRCh38.bed.gz > data/raw/rnacentral_GRCh38_coordinates.bed && rm data/raw/homo_sapiens.GRCh38.bed.gz


* The RNAcentral web interface was used to obtain the ncRNA data using the following filtering steps. The databases were downloaded to local machine and copy to the working server, uncompress and format to csv for easier manipulation.  

    * Short ncRNA: include HGNC database, exclude precursor miRNA/primary transcript, exclude rRNA, exclude lncRNA. 24/01/2024

        Query: (expert_db:"HGNC" AND TAXONOMY:"9606" AND NOT so_rna_type_name:"LncRNA" AND NOT so_rna_type_name:"Pre_miRNA" AND NOT so_rna_type_name:"RRNA" AND NOT so_rna_type_name:"MiRNA" AND length:[33 TO 19263]) AND entry_type:"Sequence"

        URL: https://rnacentral.org/search?q=expert_db:%22HGNC%22%20AND%20TAXONOMY:%229606%22%20AND%20NOT%20so_rna_type_name:%22LncRNA%22%20AND%20NOT%20so_rna_type_name:%22Pre_miRNA%22%20AND%20NOT%20so_rna_type_name:%22RRNA%22%20AND%20NOT%20so_rna_type_name:%22MiRNA%22

        gunzip data/raw/rnacentral-short-ncrna.fasta.gz
        fasta_formatter -i "rnacentral-short-ncrna.fasta" -o "rnacentral-short-ncrna-seq.csv" -t



    * Precursor miRNA: include HGNC database, include precursor miRNA/primary transcript only. 19/12/2024
        
        Query: (expert_db:"HGNC" AND TAXONOMY:"9606" AND so_rna_type_name:"Pre_miRNA" AND length:[75 TO 180]) AND entry_type:"Sequence"

        URL: https://rnacentral.org/search?q=expert_db:%22HGNC%22%20AND%20TAXONOMY:%229606%22%20AND%20so_rna_type_name:%22Pre_miRNA%22

        gunzip data/raw/rnacentral-pre-miRNA.fasta.gz
        fasta_formatter -i "rnacentral-pre-miRNA.fasta" -o "rnacentral-pre-miRNA-seq.csv" -t

    * LncRNA: include HGNC database, include lncRNA only. 28/02/2024 

         Query: expert_db:"HGNC" AND TAXONOMY:"9606" AND so_rna_type_name:"LncRNA" AND NOT so_rna_type_name:"Antisense_lncRNA" AND entry_type:"Sequence"

        URL: https://rnacentral.org/search?q=expert_db:%22HGNC%22%20AND%20TAXONOMY:%229606%22%20AND%20so_rna_type_name:%22LncRNA%22%20AND%20NOT%20so_rna_type_name:%22Antisense_lncRNA%22%20AND%20entry_type:%22Sequence%22

        gunzip data/raw/rnacentral-lncrna.fasta.gz
        fasta_formatter -i "rnacentral-lncrna.fasta" -o "rnacentral-lncrna-seq.csv" -t





#### GENOMIC REPEATS DATA
#Database for blastn
makeblastdb -in data/raw/GRCh38_p14_genomic.fna -dbtype nucl -parse_seqids -out data/raw/blastn/human_genome 

#Database for mmseqs
mkdir -p data/raw/mmseqs
mmseqs createdb data/raw/GRCh38_p14_genomic.fna data/raw/mmseqs/human_genome

#Database for mmseqs with T2T assembly
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/ 
cd data/raw/
wget https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009914755.1/download?include_annotation_type=GENOME_FASTA
# CAREFUL unzips in a folder called ncbi_dataset!
unzip download?include_annotation_type=GENOME_FASTA 
cd ../..
mv data/raw/ncbi_dataset data/raw/T2T_assembly
mkdir -p data/raw/mmseqs/T2T
mmseqs createdb "data/raw/T2T_assembly/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna" "data/raw/mmseqs/T2T/human_genome"


# Chain to map coordinates from hg38 to chm13 (i think i dont need it)
https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/
wget https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz




#### Dfam database
#data/raw
wget -P data/raw https://www.dfam.org/releases/Dfam_3.8/annotations/hg38/hg38.nrph.hits.gz
gunzip data/raw/hg38.nrph.hits.gz

# Pipe one removes the header
# Pipe two obtains Chr, Start and End columns
# Pipe three removes alternate chromosome, incomplete scaffolds and chromosomes Y/Mt


cat data/raw/hg38.nrph.hits | perl -lane 'next if(/^#|\S+_random|^chrY|^chrM|^chrUn|^chr\S+alt/);  if($F[9]>$F[10]){print "$F[0]\t$F[10]\t$F[9]"}else{print "$F[0]\t$F[9]\t$F[10]"}' | sort -k1,1V -k2,2n > data/raw/dfam-hg38-sorted.bed



##### CPC2.py
# https://github.com/gao-lab/CPC2_standalone/releases
wget -P bin/ https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz 
gzip -dc bin/v1.0.1.tar.gz | tar xf - 
mv CPC2_standalone-1.0.1 bin/
cd bin/CPC2_standalone-1.0.1
cd libs/libsvm
gzip -dc libsvm-3.18.tar.gz | tar xf -
cd libsvm-3.18
make clean && make
cd ../../../../..

# To run it
python3 bin/CPC2_standalone-1.0.1/bin/CPC2.py -i $initial_fasta -o data/cpc2-output.txt


# gnomAD database
conda install gsutil

mkdir data/raw/gnomad
cd data/raw/gnomad/
gsutil cp gs://gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr*.vcf.bgz.tbi






# ENCODE data
for file in *.bam; do
    # Check if the file exists and is a regular file
    if [ -f "$file" ]; then
        #samtools index "$file"
        total=$(samtools view -c "$file")
        echo "$file,$total" >> total-read-count.txt
    fi
done

## UNIPROT 

# Uniprot database for protein coding annotations 
10-Feb-2024
# Raw data Directory 
cd data/raw 

#   Download file
    wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/uniprot/unipAliSwissprot.bb

# Go back to working directory 
cd ../..

#   Convert bb format to bed
    ./bin/bigBedToBed data/raw/unipAliSwissprot.bb data/raw/unipAliSwissprot.bed
(chr Y and M can be removed)

#   Sort by chr first and the start coord matching order of dataset (Bedtools requirment)
    sort -k1,1V -k2,2n data/raw/unipAliSwissprot.bed > data/raw/unipAliSwissprot-sorted.bed

#   Remove extra files 
    rm -rf data/raw/unipAliSwissprot.bed
    rm -rf data/raw/unipAliSwissprot.bb


## GENCODE
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
gunzip gencode.v45.annotation.gtf.gz 

# Converting GENCODE GTF to bed file:
# Pipe two filters sequences so only known ncRNAs are included.
# Pipes three and four reformat the data into a bed file format.

cat gencode.v45.annotation.gtf | grep "Mt_rRNA\|Mt_tRNA\|miRNA\|misc_RNA\|rRNA\|scRNA\|snRNA\|snoRNA\|ribozyme\|sRNA\|scaRNA\|lncRNA" | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' > gencode-ncrna-annotation.bed

# All gencode #
cat data/raw/gencode.v45.annotation.gtf | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' > data/raw/gencode-annotation.bed

(chr Y and M can be removed)
rm -rf data/raw/gencode.v45.annotation.gtf

# GERP mammal score 
GERP scores for each sequence were extracted from the Ensembl 111-mammal GERP bigWig file using xxx (Quinlan, 2014; Yates et al., 2020).

curl -O https://ftp.ensembl.org/pub/release-111/compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw