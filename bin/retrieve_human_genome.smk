# External script for Snakefile 

###########################################################################################################################

rule download_gff:
    output:
        "data/raw/GRCh38_latest_genomic.gff"
    shell:
        "wget -O {output} https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz && gunzip {output}.gz"

rule download_fna:
    output:
        "data/raw/GRCh38_latest_genomic.fna"
    shell:
        "wget -O {output} https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz && gunzip {output}.gz"

rule reformat_fna:
    input:
        "data/raw/GRCh38_latest_genomic.fna"
    output:
        "data/raw/GRCh38_interim.csv"
    shell:
        "fasta_formatter -i {input} -o {output} -t"

rule filter_scaffolds:
    input:
        "data/raw/GRCh38_interim.csv"
    output:
        "data/raw/GRCh38.p14_genome.csv"
    shell:
        "grep 'NC_' {input} > {output} && rm {input}"
