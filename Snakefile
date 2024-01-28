# Implementing Snakemake to describe the workflow 


###########################################################################################################################

# Obtain data 
include: "retrieve_human_genome.smk" 
include: "retrieve_protein_coding_sequences.smk" 
include: non coding rna # check how implement the RNAcentral web query



rule download_rnacentral_bed:       #maybe this one on a external script 
    output:
        "data/raw/rnacentral_GRCh38_coordinates.bed"
    shell:
        """
        wget -P data/raw ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/22.0/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz &&
        gunzip -c data/raw/homo_sapiens.GRCh38.bed.gz > {output} &&
        rm data/raw/homo_sapiens.GRCh38.bed.gz
        """

#rule download_raw_sequences: 

#    output: 
#        "data/raw/rnacentral-lncrna.fasta"
#        "data/raw/rnacentral-short-ncrna.fasta"
#        "data/raw/rnacentral-pre-mirna.fasta" 


# Convert downloaded fasta files to csv for easier data retrieval 
rule fasta_to_csv:
    input: 
        lncrna="data/raw/rnacentral-lncrna.fasta",
        short_ncrna="data/raw/rnacentral-short-ncrna.fasta",
        pre_mirna="data/raw/rnacentral-pre-mirna.fasta" 
    output: 
        lncrna_csv="data/raw/rnacentral-lncrna-seq.csv", 
        short_ncrna_csv="data/raw/rnacentral-short-ncrna-seq.csv",
        pre_mirna_csv="data/raw/rnacentral-pre-mirna-seq.csv"
    conda: 
        "path to the env file" 
    shell: 
        "fasta_formatter -i {input} -o {output} -t" 


# Screening for sequences which not mitocondrial or Ychr; 75<len<3000; >3 exons 
# && extract coordinates for negative control 

rule first_screening:
    input: 
        rnacentral_coords="data/raw/rnacentral_GRCh38_coordinates.bed"
        rnacentral_lncrna="data/raw/rnacentral-lncrna-seq.csv" 
        rnacentral_short_ncrna="data/raw/rnacentral-short-ncrna-seq.csv"
        rnacentral_pre_mirna="data/raw/rnacentral-pre-mirna-seq.csv"
        genome_annotations= "data/raw/GRCh38_latest_genomic.gff"
        genome_seq="data/raw/GRCh38.p14_genome.csv"
        protein_coding_refseq="data/raw/hgnc-protein-coding-RefSeq.txt"

    output:
        protein_exon_two_csv="data/protein-exon2-dataset.csv"
        protein_exon_three_csv="data/protein-exon3-dataset.csv"
        lncrna_exon_two_csv="data/functional-lncrna-exon2-dataset.csv"
        lncrna_exon_three_csv="data/functional-lncrna-exon3-dataset.csv"
        short_ncrna_csv="data/functional-short-ncrna-dataset.csv"
        pre_mirna_csv="data/functional-pre-mirna-dataset.csv"
        negative_control_csv="data/coords-for-negative-control.csv"
    conda:
    
    script:
        bin/initial-dataset.sh


# Generate functional FASTA file from CSV file
names=["protein-exon2", "protein-exon3", "functional-lncrna-exon2", "functional-lncrna-exon3", "functional-short-ncrna", "functional-pre-mirna"]

rule all:
    input:
        expand("data/{name}-seq.fa", name=names)

rule csv_to_fasta:
    input:
        csv="data/{name}-dataset.csv"
    output:
        fasta="data/{name}-seq.fa"
    shell:                                      # should this be in python? 
        """
        while IFS=',' read -r id _ _ _ _ seq
        do
            echo -e ">$id\n$seq" >> {output.fasta}
        done < {input.csv}
        """