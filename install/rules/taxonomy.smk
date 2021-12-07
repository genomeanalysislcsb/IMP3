##################################################
# Kraken2 databases
##################################################
rule kraken2db_download:
    output:
        download=temp(directory("downloads/files/{kdb}")),
    log:
        "{kdb}/download.log"
    wildcard_constraints:
        kdb="|".join(config["kraken2"].keys())
    params:
        url=smk_get_kraken2db_url
    conda:
        os.path.join(ENV_DIR, "utils.yaml")
    message:
        "Kraken2 database download: {params.url}"
    script:
        os.path.join(SRC_DIR, "download_and_decompress.py")

# NOTE: folder structure of the downloaded and decompressed archive is NOT preserved
rule kraken2db_files:
    input:
        "downloads/files/{kdb}" # folder (downloaded and decompressed archive)
    output:
        files=expand("{{kdb}}/{fname}", fname=get_kraken2db_file_names()),
        flist=temp("{kdb}/files.txt"),
    log:
        "{kdb}/{kdb}.log"
    wildcard_constraints:
        kdb="|".join(config["kraken2"].keys())
    message:
        "Kraken database files"
    shell:
        "find {input}/ -type f ! -name '.snakemake_timestamp' ! -name 'imp3_install_*' | sed 's@{input}/@@' | sort > {output.flist} && "
        "rsync -rtP --no-R --log-file={log} --files-from={output.flist} {input} $(dirname {output.flist})/ 2>> {log}"

##################################################
# GTDB-tk
##################################################

rule gtdbtk_download:
    output:
        result=directory("gtdbtk"),
        download=temp(directory("downloads/files/gtdbtk")),
    log:
        "gtdbtk/download.log"
    params:
        url=config["gtdbtk"]["url"]
    conda:
        os.path.join(ENV_DIR, "utils.yaml")
    message:
        "GTDB-tk data download: {params.url}"
    script:
        os.path.join(SRC_DIR, "download_and_decompress.py")
