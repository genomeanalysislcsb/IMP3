##################################################
# Adapters
##################################################

# Adapters from Trimmomatic repository
rule trimmomatic_adapters:
    input:
        "downloads/git/trimmomatic"
    output:
        os.path.join("adapters", "trimmomatic.done")
    log:
        os.path.join(INSTALL_DIR, "logs/trimmomatic_setup.log")
    message:
        "Adapters: Trimmomatic adapters"
    shell:
        "rsync -rtP --log-file={log} {input}/adapters/ $(dirname {output})/ && touch {output}"

##################################################
# rRNA
##################################################

# rRNA references from SortMeRNA repository
rule sortmerna_rrna:
    input:
        "downloads/git/sortmerna"
    output:
        os.path.join("sortmerna", "sortmerna.done")
    log:
        os.path.join(INSTALL_DIR, "logs/sortmerna_setup.log")
    message:
        "Filtering: SortMeRNA files"
    shell:
        "rsync -rtP --log-file={log} {input}/data/rRNA_databases/*.fasta $(dirname {output})/ && touch {output}"

##################################################
# Genomes
##################################################
rule genome_download:
    output:
        result="filtering/{genome}.fa",
        download=temp(directory("downloads/filtering/{genome}")),
    log:
        os.path.join(INSTALL_DIR, "logs/{genome}_download.log")
    wildcard_constraints:
        genome="|".join(config["genomes"].keys())
    params:
        url=smk_get_genome_url
    conda:
        os.path.join(ENV_DIR, "utils.yaml")
    message:
        "Filtering: download genome: {params.url}"
    script:
        os.path.join(SRC_DIR, "download_and_decompress.py")