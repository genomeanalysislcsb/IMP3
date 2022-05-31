##################################################
# Kraken2 databases
##################################################
# NOTE: Expected dtabase files: hash.k2d, opts.k2d, taxo.k2d
#       https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases
# NOTE: The folder structure of the downloaded and decompressed archive will NOT be preserved

rule kraken2db_download:
    output:
        download=temp(directory("downloads/files/{kdb}")),
    log:
        os.path.join(INSTALL_DIR, "logs/{kdb}_download.log")
    wildcard_constraints:
        kdb="|".join(config["kraken2"].keys())
    params:
        url=smk_get_kraken2db_url
    threads: 4
    conda:
        os.path.join(ENV_DIR, "utils.yaml")
    message:
        "Kraken2 database download: {params.url}"
    script:
        os.path.join(SRC_DIR, "download_and_decompress.py")

rule kraken2db_files:
    input:
        "downloads/files/{kdb}" # folder (downloaded and decompressed archive)
    output:
        files=expand("{{kdb}}/{fname}", fname=get_kraken2db_file_names()),
        flist=temp("{kdb}/files.txt"),
    log:
        os.path.join(INSTALL_DIR, "logs/{kdb}_setup.log")
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
        download=temp(directory("downloads/files/gtdbtk"))
    log:
        os.path.join(INSTALL_DIR, "logs/gtdbtk_download.log")
    params:
        url=config["gtdbtk"]["url"]
    threads:
        4
    conda:
        os.path.join(ENV_DIR, "utils.yaml")
    message:
        "GTDB-tk data download: {params.url}"
    script:
        os.path.join(SRC_DIR, "download_and_decompress.py")