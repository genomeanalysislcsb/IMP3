##################################################
# GTDB-tk
##################################################

rule gtdbtk_download:
    output:
        download = directory("gtdbtk")
    log:
        os.path.join(INSTALL_DIR, "logs/motus_download.log")
    threads:
        4
    params:
        url = config["gtdbtk"]["url"]
    conda:
        os.path.join(ENV_DIR, "utils.yaml")
    message:
        "GTDB-tk data download: {params.url}"
    script:
        os.path.join(SRC_DIR, "download_and_decompress.py")
