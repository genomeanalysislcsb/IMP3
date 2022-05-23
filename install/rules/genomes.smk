##################################################
# Genomes
##################################################
rule genome_download:
    output:
        result="filtering/{genome}.fa",
        download=temp(directory("downloads/filtering/{genome}")),
    log:
        "filtering/{genome}.log"
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