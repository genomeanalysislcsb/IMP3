##################################################
# GTDB-tk
##################################################

rule gtdbtk_download:
    output:
        db_dir=directory("gtdbtk")
    log:
        os.path.join(INSTALL_DIR, "logs/motus_download.log")
    threads:
        4
    conda:
        os.path.join(ENV_DIR, "utils.yaml")
    message:
        "GTDB-tk data download: {params.url}"
    script:
        os.path.join(SRC_DIR, "download_and_decompress.py")
