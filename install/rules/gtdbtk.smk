##################################################
# GTDB-tk
##################################################

rule gtdbtk_download:
    output:
        db_dir = directory("gtdbtk"),
        final_output = directory("gtdbtk/release207_v2")
    log:
        os.path.join(INSTALL_DIR, "logs/gtdbtk_download.log")
    params:
        url = config["gtdbtk"]["url"],
        fname = config["gtdbtk"]["url"].split('/')[-1]
    threads:
        4
    conda:
        os.path.join(ENV_DIR, "utils.yaml")
    message:
        "GTDB-tk data download: {params.url}"
    shell:
        """
        wget -P {output.db_dir} {params.url} >> {log} 2>&1 && \
        pigz -d -c {output.db_dir}/{params.fname} | tar -x -f - -C {output.db_dir} >> {log} 2>&1
        """
