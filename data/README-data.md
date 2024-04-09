# DATA - Features of Functional Human Genes 

* [Databases for Sequences Datasets](#databases-for-sequences-datasets)
    * [Human Genome](#human-genome)
    * [Protein-coding Genes](#protein-coding-genes)
    * [Non-coding Genes](#non-coding-genes)
        * [lncRNA](#lncrna)
        * [short-ncRNA](#short-ncrna)
        * [pre-miRNA](#precursor-mirna)
    * [Negative Control](#negative-control)
        * [Conservation Scores](#conservation-scores)
        * [Transcriptome Expression](#transcriptome-expression)
        * [Genomic Repeats](#genomic-repeats)
        * [Specific RNA and Protein-Coding ](#specific-rna-and-protein-coding)
        * [Population Variation](#population-variation)

* [References](#references)

# Raw Data: 

The following databases and files should be stored on the data/raw directory

    cd data/raw

# Databases for Sequences Datasets
## Human Genome: 

The human genome is obtain from xxxx 

* #### Human genome annotations (Last access: 15/12/2023)

    Download and unzip GFF file

        wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz   
    
        gunzip GRCh38_latest_genomic.gff.gz

        
* ####  Human genome sequence (Last access: 15/12/2023)
    
    Download and unzip FNA file 

        wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz   && 
    
        gunzip GRCh38_latest_genomic.fna.gz

    Format file to a tab delimited CSV file using the FASTA formatter from [FASTX-Toolkit version 0.0.13](http://hannonlab.cshl.edu/fastx_toolkit/) removing scaffolds and alternative chromosomes, to allow for easier manipulation and parsing

        fasta_formatter -i GRCh38_p14_genomic.fna -o GRCh38_interim.csv -t 

        grep "NC_" GRCh38_interim.csv > GRCh38.p14_genome.csv
        
        rm GRCh38_interim.csv

 
## Protein-coding Genes:   

The IDs for the protein-coding genes are obtain from the Hugo Gene Nomenclature Committee (HGNC) (Braschi et al., 2019)

* #### Generate a file with RefSeq IDs of protein-coding genes: 

    Download text file for all protein-coding genes from HGNC  Last access: 15/12/2023
    
        wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt

    Process file to remove sequences encoded in the Y and Mt chromosomes or are no longer approved genes (pipes 1&2), and extract columns associated with RefSeq IDs, and removes sequences associated with more than one ID (pipes 3-5)
   
        grep -v "mitochondrially encoded\|Entry Withdrawn\|Yq" gene_with_protein_product.txt | grep "Approved" | cut -f 24 | grep -v "|" | grep . > hgnc-protein-coding-RefSeq.txt


## Non-coding Genes:

The Non-coding genes coordinates are obtain from RNAcentral (The RNAcentral Consortium, 2019)

* #### Chromosome coordinates for all non-coding RNA:

    Download the Homo sapiens GRCh38 bed file from the RNAcentral FTP directory, which contains the chromosome coordinates for all human ncRNAs; and filter out patches from RNAcentral coordinates database keeping only coordinates corresponding to chr1-22 and X (Last access: 19/12/2024)

        wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/22.0/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz

        gunzip -c homo_sapiens.GRCh38.bed.gz | sort -k1,1V -k2,2n | grep -E '^chr(1?[0-9]|2[0-2]|X)\b' | awk -F'\t' '!seen[$1,$2,$3]++' > rnacentral-GRCh38-coords.bed && rm homo_sapiens.GRCh38.bed.gz


* #### Non-coding RNA annotations 

    Use the [RNAcentral web interface](https://rnacentral.org/) to obtain the ncRNA annotations using the following filtering steps. Download the databases to local machine and copy them to the working directory, uncompress and format to csv for easier manipulation.  

    * #### lncRNA 
    Include HGNC database, include lncRNA only (Last access: 28/02/2024) 

        Query: expert_db:"HGNC" AND TAXONOMY:"9606" AND so_rna_type_name:"LncRNA" AND NOT so_rna_type_name:"Antisense_lncRNA" AND entry_type:"Sequence"  --> 2,898 sequences

        URL: https://rnacentral.org/search?q=expert_db:%22HGNC%22%20AND%20TAXONOMY:%229606%22%20AND%20so_rna_type_name:%22LncRNA%22%20AND%20NOT%20so_rna_type_name:%22Antisense_lncRNA%22%20AND%20entry_type:%22Sequence%22

        gunzip rnacentral-lncrna.fasta.gz

        fasta_formatter -i "rnacentral-lncrna.fasta" -o "rnacentral-lncrna-seq.csv" -t

    * #### Short ncRNA
    Include HGNC database, exclude precursor miRNA/primary transcript, exclude rRNA, exclude lncRNA (Last access: 24/01/2024)

        Query: (expert_db:"HGNC" AND TAXONOMY:"9606" AND NOT so_rna_type_name:"LncRNA" AND NOT so_rna_type_name:"Pre_miRNA" AND NOT so_rna_type_name:"RRNA" AND NOT so_rna_type_name:"MiRNA" AND length:[33 TO 19263]) AND entry_type:"Sequence"

        URL: https://rnacentral.org/search?q=expert_db:%22HGNC%22%20AND%20TAXONOMY:%229606%22%20AND%20NOT%20so_rna_type_name:%22LncRNA%22%20AND%20NOT%20so_rna_type_name:%22Pre_miRNA%22%20AND%20NOT%20so_rna_type_name:%22RRNA%22%20AND%20NOT%20so_rna_type_name:%22MiRNA%22

        gunzip rnacentral-short-ncrna.fasta.gz

        fasta_formatter -i "rnacentral-short-ncrna.fasta" -o "rnacentral-short-ncrna-seq.csv" -t

    * #### Precursor miRNA
    Include HGNC database, include precursor miRNA/primary transcript only (Last access: 19/12/2024)
        
        Query: (expert_db:"HGNC" AND TAXONOMY:"9606" AND so_rna_type_name:"Pre_miRNA" AND length:[75 TO 180]) AND entry_type:"Sequence"

        URL: https://rnacentral.org/search?q=expert_db:%22HGNC%22%20AND%20TAXONOMY:%229606%22%20AND%20so_rna_type_name:%22Pre_miRNA%22

        gunzip rnacentral-pre-miRNA.fasta.gz

        fasta_formatter -i "rnacentral-pre-miRNA.fasta" -o "rnacentral-pre-miRNA-seq.csv" -t


## Negative control
The chromosome coordinates and sequences for the final negative control protein-coding, short ncRNA and lncRNA genes are available in [data](data/) (ADD LINK TO GITHUB!!!)

The coordinates for the negative control datasets are extracted from regions of the genome lacking known annotated genes. Both GENCODE and RNAcentral databases are used to determinate this genes complement regions. 

#### GENCODE database

Download GENCODE annotations of known genes (Last access: 13/11/2023 - version 45):

    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz

    gunzip gencode.v45.annotation.gtf.gz 

Converting GENCODE GTF to bed file, including only genes and removing chrM/Y:

    cat gencode.v45.annotation.gtf | awk 'OFS="\t" {if ($3=="gene" && $1 != "chrM" && $1 != "chrY") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' > gencode-annotation.bed

Remove excess file:

    rm -rf gencode.v45.annotation.gtf


#### Genes Complement 

Obtain genes complement using bedtools complement. The length of the human chromosomes was used as genome reference. 

    gencode_bed="data/raw/gencode-annotation.bed"
    rnacentral_bed="data/raw/rnacentral-GRCh38-coords.bed"
    chromo_len="data/raw/chromosomes-len.bed"

    bedtools complement -i "$gencode_bed" -g "$chromo_len" > data/raw/gencode-complement
    bedtools complement -i "$rnacentral_bed" -g "$chromo_len" > data/raw/rnacentral-complement

Merge both files, sort and remove duplicated entries

    cat data/raw/gencode-complement data/raw/rnacentral-complement | sort -k1,1V -k2,2n | uniq > data/raw/genes-complement.bed && rm data/raw/rnacentral-complement && rm data/raw/gencode-complement


# Databases for Features 

## Conservation Scores
#### phyloP 241-mammal 
Download the bigWig file from the Zoonomia project resouces (reference) (Last access: 11/12/2023)

    wget https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Homo_sapiens/241-mammalian-2020v2.bigWig 

### phastCons 
--> probably will be removed it from project 

#### GERP mammal score 
Download the bigWig file from Ensembl 111-mammal resource (Quinlan, 2014; Yates et al., 2020 ????) (Last access: 30/12/2023)

    curl -O https://ftp.ensembl.org/pub/release-111/compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw


## Transcriptome expression
The links for the ENCODE total and small RNA-Seq data (BAM files) for [primary-cell](data/encode-primary-cell-rnaseq.txt) and [tissue](data/encode-tissue-rnaseq.txt) samples are available in the [data](data/) folder (ENCODE Project Consortium, 2012), with all the RNA-Seq data kept in sample specific folders. 
The BAM files need to be indexed using ```samtools index``` and the total reads per file calculated into a text file in the same folder (Li et al., 2009).

Download ENCODE RNA-Seq datasets (must be done within a specified folder for each sample type)

    mkdir ENCODE-tissue && mkdir ENCODE-primary-cell
    
    cd ENCODE-tissue/

    xargs -L 1 curl -O -J -L < ../encode-tissue-rnaseq.txt         # Done within ENCODE-tissue/

If your downloads fail, delete empty files and restart downloads (check your rsync/backups aren't duplicating these!!!):

    ls -s *bam | perl -lane 'print "rm $F[1]" if($F[0]==0)' | sh && cat ../encode-tissue-rnaseq.txt | perl -lane 'if(/download\/(ENC\S+.bam)/){ if(-s $1){print "#$_"}else{print "$_"}  }' > blah && mv blah ../encode-tissue-rnaseq.txt

    cd ../ENCODE-primary-cell/
    xargs -L 1 curl -O -J -L < ../encode-primary-cell-rnaseq.txt   # Done within ENCODE-primary-cell/

If your downloads fail, delete empty files and restart downloads:

    ls -s *bam | perl -lane 'print "rm $F[1]" if($F[0]==0)' | sh && cat ../encode-primary-cell-rnaseq.txt | perl -lane 'if(/download\/(ENC\S+.bam)/){ if(-s $1){print "#$_"}else{print "$_"}  }' > blah && mv blah ../encode-primary-cell-rnaseq.txt

    nprocs=12
    for url in $(cat to-get.txt); 
        do 
        curl -O -J -L $url & 
        pwait "$nprocs"
    done

Index files and calculate total reads per RNA-Seq dataset (done within each specified sample type folder [ENCODE-tissue/] & [ENCODE-primary-cell/])

    for file in *.bam; do
    # Check if the file exists and is a regular file
    
        if [ -f "$file" ]; then
            #samtools index "$file"
            total=$(samtools view -c "$file")
            echo "$file,$total" >> total-read-count.txt
        fi
    done

PPG: for indexing partially downloaded/indexed runs: 

    ls *bai | perl -lane 's/.bai$//g; print' | tr "\n" "|"
    ls *bam | egrep -v 'ENCFF003XND.bam|ENCFF015GPO.bam|ENCFF018MLY.bam|ENCFF022KDS.bam|ENCFF025WMR.bam|ENCFF040ULF.bam|ENCFF044SZR.bam|ENCFF053EDJ.bam|ENCFF054ZIN.bam|ENCFF063VHB.bam|ENCFF072WPN.bam|ENCFF085BXZ.bam|ENCFF096FNU.bam|ENCFF124IFE.bam|ENCFF136HDU.bam|ENCFF139VRN.bam|ENCFF146SZI.bam|ENCFF151OFT.bam|ENCFF151ZTF.bam|ENCFF157HFY.bam|ENCFF164UED.bam|ENCFF165SZN.bam|ENCFF166EQZ.bam|ENCFF176XWZ.bam|ENCFF199BEE.bam|ENCFF227RIN.bam|ENCFF263ZZQ.bam|ENCFF265AFD.bam|ENCFF265UEF.bam|ENCFF272UGL.bam|ENCFF274DJM.bam|ENCFF276QUS.bam' | perl -lane 'print "samtools index $F[0] && total=\$( samtools view $F[0] | wc -l ) && echo \$file,\$total >> total-read-count.txt" ' | sh

    ls *bam | egrep -v 'ENCFF017ZWH.bam|ENCFF019KCS.bam|ENCFF023PAK.bam|ENCFF029BEL.bam|ENCFF058KEP.bam|ENCFF071GAO.bam|ENCFF090FOY.bam|ENCFF116PFD.bam|ENCFF119POR.bam|ENCFF128NQX.bam|ENCFF159QAX.bam|ENCFF164JEV.bam|ENCFF166RKD.bam|ENCFF169KJY.bam|ENCFF184YUO.bam|ENCFF196YLT.bam|ENCFF202ETQ.bam|ENCFF220CVY.bam|ENCFF224PNI.bam|ENCFF333FFK.bam|ENCFF335UAF.bam|ENCFF376UOL.bam|ENCFF403LXD.bam|ENCFF406DEH.bam|ENCFF407VTT.bam|ENCFF413CTV.bam|ENCFF422EVZ.bam|ENCFF444EKY.bam|ENCFF460PJA.bam|ENCFF491VZY.bam|ENCFF514QIL.bam|ENCFF542GSD.bam|ENCFF549UBP.bam|ENCFF559BDS.bam|ENCFF575PML.bam|ENCFF584FBF.bam|ENCFF590PJR.bam|ENCFF592AWU.bam|ENCFF595LQV.bam|ENCFF596OPC.bam|ENCFF604QYU.bam|ENCFF681WQP.bam|ENCFF687MKT.bam|ENCFF701CAC.bam|ENCFF730CJE.bam|ENCFF741AOJ.bam|ENCFF751DPV.bam|ENCFF756NWQ.bam|ENCFF778XJG.bam|ENCFF804YQY.bam|ENCFF806MFZ.bam|ENCFF811MSS.bam|ENCFF816ZAP.bam|ENCFF836WKP.bam|ENCFF840PMN.bam|ENCFF841GMR.bam|ENCFF841MML.bam|ENCFF867ROO.bam|ENCFF883NIZ.bam|ENCFF969GXB.bam|ENCFF979LYG.bam' | perl -lane 'print "samtools index $F[0] && total=\$( samtools view $F[0] | wc -l ) && echo \$file,\$total >> total-read-count.txt" '


## Genomic Repeats

* #### Copy number: 
Generate a local database from the previously downloaded human genome (add link refering above) sequence required to run mmseqs tool (Altschul et al., 1990; O’Leary et al., 2016???) 

    mkdir -p mmseqs
    mmseqs createdb GRCh38_p14_genomic.fna mmseqs/human_genome

* #### Distance to the nearest repetitive element: 
Generate a file with Dfam non-redundent hg38 repetitive DNA elements.

Download Dfam database: 

    wget -P https://www.dfam.org/releases/Dfam_3.8/annotations/hg38/hg38.nrph.hits.gz
    
    gunzip hg38.nrph.hits.gz

Format file to removes the header (pipe 1), obtain Chr, Start and End columns (pipe 2) and remove alternate chromosome, incomplete scaffolds and chromosomes Y/Mt (pipe 3):

    cat hg38.nrph.hits | perl -lane 'next if(/^#|\S+_random|^chrY|^chrM|^chrUn|^chr\S+alt/);  if($F[9]>$F[10]){print "$F[0]\t$F[10]\t$F[9]"}else{print "$F[0]\t$F[9]\t$F[10]"}' | sort -k1,1 -k2,2n -k3,3n > dfam-hg38-sorted.bed


## Population variation 

Download SNP data from the Genome Aggregation Database (gnomAD) for each chromosome. For this, use gsutil (conda install gsutil) xxx 

    mkdir gnomad
    
    cd gnomad/
    
    gsutil cp gs://gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr*.vcf.bgz.tbi



## Go back to Working directory 

    cd ../..






