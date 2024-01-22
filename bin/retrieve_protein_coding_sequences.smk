
# External script for Snakefile

###########################################################################################################################

rule download_hgnc_data:
    output:
        "data/raw/gene_with_protein_product.txt"
    shell:
        "wget -O {output} ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt"

rule filter_protein_coding_genes:
    input:
        "data/raw/gene_with_protein_product.txt"
    output:
        "data/raw/hgnc-protein-coding-RefSeq.txt"
    shell:
        """
        grep -v "mitochondrially encoded\|Entry Withdrawn\|Yq" {input} | \
        grep "Approved" | \
        cut -f 24 | \
        grep -v "|" | \
        grep . > {output}
        """

