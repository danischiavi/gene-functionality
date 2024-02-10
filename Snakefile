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



########################

# TO EXTRACT THE FEATURES 

## GC content 

rule all_GC:
    input:
        expand("data/{sample}-GC.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"])

rule run_intrinsic_seq_features:
    input:
        "data/{sample}-dataset.csv"
    output:
        "data/{sample}-GC.csv"
    shell:
        "bin/intrinsic-seq-features.sh {input} > {output}"


## Conservation 

rule all_conservation:
    input:
        expand("data/conservation/{sample}-conservation.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"])

rule run_seq_conservation_features:
    input:
        "data/{sample}-dataset.csv"
    output:
        "data/conservation/{sample}-conservation.csv"
    shell:
        "bin/seq-conservation-features.sh {input} > {output}"


## Transcriptome

rule all_transcriptome:
    input:
        expand("data/transcriptome/{sample}-transcriptome.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"])

rule run_transcriptome:
    input:
        "data/{sample}-dataset.csv"
    output:
        "data/transcriptome/{sample}-transcriptome.csv"
    shell:
        "bin/transcriptome-expression-features.sh {input} > {output}"


## Genomic repeats 

rule all_genomic_repeat:
    input:
        expand("data/repeats/{sample}-copy-number.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"]),
        expand("data/repeats/{sample}-dfam-distance.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"])

rule run_genomic_repeat:
    input:
        initial_data="data/{sample}-dataset.csv",
        initial_fasta="data/{sample}-seq.fa"
    output:
        "data/repeats/{sample}-copy-number.csv",
        "data/repeats/{sample}-dfam-distance.csv"
    shell:
        "bin/genomic-repeat-associated-features.sh {input} > {output}"


## Population variation

rule all_population_variation:
    input:
        expand("data/population/{sample}-1kGP-variation.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"]),
        expand("data/population/{sample}-gnomAD-variation.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"])

rule run_populatio_variation:
    input:
        initial_data="data/{sample}-dataset.csv",
        initial_fasta="data/{sample}-seq.fa"
    output:
        "data/population/{sample}-1kGP-variation.csv",
        "data/population/{sample}-gnomAD-variation.csv"
    shell:
        "bin/population-variation-features.sh {input} > {output}"


## protein and RNA specific 

rule all_protein_and_rna_specific_features:
    input:
        expand("data/specific/{sample}-interaction.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"]),
        expand("data/specific/{sample}-structure.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"]),
        expand("data/specific/{sample}-coding-potential.csv", sample=["protein-exon2", "protein-exon3", "functional-lncrna-exon1", "functional-lncrna-exon2", "functional-short-ncrna"])


rule run_protein_and_rna_specific_features:
    input:
        initial_data="data/{sample}-dataset.csv",
        initial_fasta="data/{sample}-seq.fa"
    output:
        "data/specific/{sample}-interaction.csv",
        "data/specific/{sample}-structure.csv",
        "data/specific/{sample}-coding-potential.csv"
    shell:
        "bin/protein-and-rna-specific-features.sh {input} > {output}"