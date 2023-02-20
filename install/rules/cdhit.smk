##################################################
# CD_HIT
##################################################

rule cdhit_compile:
    input:
        "mantis/mantis.config"
    output:
        cdhit_dir = directory(os.path.join(SUBMOD_DIR, "cdhit")),
        cdhit_binary = os.path.join(SUBMOD_DIR, "cdhit", "cd-hit-est")
    log:
        os.path.join(INSTALL_DIR, "logs/cdhit_setup.log")
    params:
        url = config["cdhit"]["url"],
        max_seq_len = config["cdhit"]["max_seq_len"],
        fname = config["cdhit"]["url"].split("/")[-1]
    threads:
        workflow.cores
    conda:
        os.path.join(ENV_DIR, "utils.yaml")
    message:
        "CD_HIT: Compiling CD_HIT with max seq lenght {params.max_seq_len}"
    shell:
        """
        wget -nv {params.url} \
             -P {output.cdhit_dir} >> {log} 2>&1 && \
        tar -xvf {output.cdhit_dir}/{params.fname} \
            -C {output.cdhit_dir} \
            --strip-components=1 >> {log} 2>&1 && \
        rm -rf {output.cdhit_dir}/{params.fname} && \
        cd {output.cdhit_dir} && \
        make MAX_SEQ={params.max_seq_len} >> {log} 2>&1
        """
