##################################################
# GTDB-tk
##################################################

rule gtdbtk_download:
    output:
        result=directory("gtdbtk"),
        download=temp(directory("downloads/files/gtdbtk")),
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
