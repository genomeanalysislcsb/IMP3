include:
    "../../Multi/Catalogue/cluster.smk"
if config['catalogue']['kraken']['do']:
    include:
        "../../Multi/Catalogue/kraken_genes.smk"
if config['catalogue']['mapping']['do']:
    include:
        "../../Multi/Catalogue/mapping.smk"
include:
    "../../Multi/Catalogue/table.smk"


# master command
rule all_cat:
    input:
        'Catalogue/catalogue.faa',
        'Catalogue/aa_catalogue.faa',
        'Catalogue/catalogue.fasta',
        'Catalogue/catalogue.tsv'
    output:
        touch('status/catalogue.done')


