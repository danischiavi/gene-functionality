# Features of Functional Human Genes 

 This is the repository for the "Features of Functional Human Genes" paper pipeline, data and results.

## Simple access instructions

To access the project, create a local project directory. For example:

    mkdir gene-functionality 

Clone this repository to the project directory.

    git clone https://github.com/danischiavi/gene-functionality gene-functionality

This will download the project pipline and curated paper data and results. The data and results corresponding to the paper are located in the XXX directory.

If you want to use any of the scripts, you will need to perform further installation steps.


## Dependencies

Most dependencies have been set in the conda environemnt YML file.

Create a new conda environment using the YML file loated in the scripts directory.

    conda env create -f condaEnv.yml

Some software was unavailable for adding to the conda environment and will need to be installed/downloaded manually into the **./bin** project directory 

 * #### UCSC binary utilities 

This directory contains Genome Browser and Blat application binaries built for standalone command-line use on various supported Linux and UNIX platforms
https://hgdownload.soe.ucsc.edu/admin/exe/

* bigBedToBed

        wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed 

 * bigWigSummary

        wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary 

 

## Pipeline Purpose

This pipeline has been created to evaluate XXX. For the research paper that is relevant to this pipeline, please see: XXX


## How to use this pipeline

When running the pipeline, always run from the **project directory**, unless otherwise specified.

1) Download required Databases as descripted in [data/README-data.md file](LINK)  

2) Generate functional and negative control databases following the [Snakefile-sequences file](LINK).For this, run the following command on the Working directory of the project replacing "cores" for the number of cores you want to use  

        snakemake -s Snakefile-sequences -j "cores" 

3) Obtain features for each dataset following [Snakefile-features file](LINK). For this, run the following command on the Working directory of the project replacing "cores" for the number of cores you want to use  

        snakemake -s Snakefile-features -j "cores"

4) Analyze results with [scripts/spearman-function-feature.R](link) and [scripts/random-forest-analysis.R](link) for the main analysis. For complementary analysis see the [complementary directory](link)


## Description of scripts

* #### A0-protein-coding-sequences.sh

Input:  RefSeq IDs from HGNC for coding potential RNA sequences, human genome annotations and sequence
Output: functional coding protein RNA sequences datasets (exon 2&3) and coordinates for their negative control

Generates the functional protein-coding RNA dataset for exon 2 and 3 with 1000 sequences, and the coordinates to generate their negative control dataset. For this, we randomly select protein-coding genes IDs from the HGNC database (linkxx) and extract the corresponding coordinates and sequences from the GRCh38 genome (annotations and fasta file linkxxx ncbi). We exclude sequences by chromosome Y and mitrochondria, exon number, sequence length and unknown nucleotides (see maniscript for details) 

For the negative coordinates, we extract the start coordinate of the first exon and end coordinate of last exon, and the length of exon 2 and 3. This files is used [downstream](#A3-negative-control.sh) as input to generate the negative control dataset


* #### A1-lncRNA-sequences.sh

Input: lncRNA sequences and non-coding RNA annotations from RNAcentral 
Output: functional lncRNA sequences datasets (exon 1&2) and coordinates for their negative control

Generates the functional lncRNA dataset for exon 1 and 2 with 1000 sequences, and the coordinates to generate their negative control dataset. For this, we randomly select lncRNA genes IDs from the RNAcentral database (linkxx) and extract the corresponding coordinates and sequences from the RNAcentral as well (annotations and fasta file linkxxx). We filter out sequences and generate the coordinates for the negative control as described for the functional protein-coding dataset.   

* #### A2-short-ncRNA-sequences.sh

Input: short-ncRNA sequences and non-coding RNA annotations from RNAcentral 
Output: functional short-ncRNA sequences datasets (exon 1&2) and coordinates for their negative control

Generates the functional short-non-coding RNA dataset with 1000 sequences, and the coordinates to generate its negative control dataset as described for the functional lncRNA dataset. 

* #### A3-negative-control.sh 

Input: gene coordinates and genome genes complement (genome regions lacking of genes - for details see)
Output: negative control dataset for each functional dataset 

Generates the negative control dataset of each functional dataset. For this, we input each file with the coordinates for the negative control generated as described above (Ex: xxxx). From each exon on this file, we generate 10 other coordinates (5 downstream and 5 upstream) considering the length of the sequence and a fixed distance from the corresponding gene.   
We search for the closest genomic region of the genes complement to each of these generated coordinates. These new coordinates are used to extract the corresponding sequences from the human genome and generate the final negative control dataset.  


* #### B0-intrinsic-seq-features.sh

Calculates the percentage of G+C content of each sequence, counting the number of each nucleotide. Also determinates the low complexity sequences density. For this, we identify low sequence regions within each sequences and count the corresponding number of nucleotides to calculate the density.  

* #### B1-seq-conservation-features.sh

Extracts the maximum PhyloP conservation score for each sequence from Zoonomia resources. 

* #### B2-transcriptome-expression-features.sh

Calculates the RPKM (average number of reads per kilobase of transcripts per million mapped reads) for each sequence. For this, the reads and total reads of small RNA-Seq datasetsfrom 71 human tissue and 140 human primary-cell samples are obtained from ENCODE. 

* #### B3-genomic-repeats-associated-features.sh

Calculates the genomic copy number of each sequence on the human genome version GRCh38. Also, calculates the distance to the closest non-overlapping transposable DNA elements in the human genome obtained from Dfam vxxx.

* #### B4-protein-and-rna-specific-features.sh

Calculates coding potential using RNAcode based on a multiple sequence alignment 241way from Zoonomia resoures.  

Calculates maximum covariance score to analyse the structure. Also, calculates the Minumum Free Energy.  

Calculates the average interaction free energy between the sequences and a curated database of 34 human ncRNA sequences from RNAcentral v15

* #### B5-population-variation-features.sh

Calculates Minor Allele Frequency (MAF) and the Single Nucleotide Polymorphisms (SNPs) density. For this, the SNPs count for each sequence is extracted from gnomAD v3 database. 

* #### C0-spearman-function-feature.R

Calculates a Spearman correlation between each feature and functionality. The correlation matrix for each dataset is generated with the confidence intervals calculated using a bootstrap approach (N=1,000).

* #### C1-random-forest-analysis.R

Generates a Random forests with 1,000 trees, using 70% training data and 30% test data. A bootstrap approach is used to simulate a further 100 classification models to assess performance variation due to sampling biases. 

#### Others 

* B4.1-maf-to-stk.sh

    Converts multiple sequence alignment MAF files to Stockholm. This is used withing the protein-and-rna-specific-features.sh script. 


