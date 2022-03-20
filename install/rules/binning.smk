##################################################
# Checkm setup for MetaWrap
##################################################

rule setup_checkm:
    output:
        checkm_db=directory("checkm"),
        setup_check="checkm/checkm_db.ready"
    conda:
        os.path.join(MAIN_ENV_DIR, "IMP_metawrap.yaml")
    log:
        os.path.join(INSTALL_DIR, "logs/MetaWrap_CheckM_setup.log")
    message:
        "MetaWrap: Setting up CheckM."
    shell:
        """
        echo "Downloading CheckM database" >> {log} 2>&1
        wget -nv https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -P {output.checkm_db} >> {log} 2>&1 && \
        echo "Decompressing CheckM database"  >> {log} 2>&1&& \
        tar -xvf {output.checkm_db}/checkm_data_2015_01_16.tar.gz -C {output.checkm_db} >> {log} 2>&1 && \
        rm -rf {output.checkm_db}/checkm_data_2015_01_16.tar.gz &&\
        echo "Setting up CheckM" >> {log} 2>&1 &&\
        printf "%s\n" {output.checkm_db} | checkm >> {log} 2>&1 && touch {output.setup_check} && \
        echo "Setup done." >> {log} 2>&1
        """

##################################################
# binny
##################################################

rule binny_setup:
    input:
        "mantis/setup.done"
    output:
        binny_dir=directory(os.path.join(SUBMOD_DIR, "binny")),
        out=os.path.join(SUBMOD_DIR, "binny/setup.done")
    log:
        os.path.join(INSTALL_DIR, "logs/binny_setup.log")
    params:
        prokka_env = os.path.join(MAIN_ENV_DIR,"IMP_annotation.yaml"),
        mantis_env=os.path.join(MAIN_ENV_DIR, "IMP_mantis.yaml"),
        snakemake_env=os.path.join(CONDA_DIR, "snakemake_env")
    threads:
        4
    message:
        "binny: Setting up binny"
    shell:
        """
        echo "Downloading binny." >> {log} 2>&1
        git clone -b binny_beta_v0.50_multi_bam_v031 --single-branch https://github.com/ohickl/binny.git \
        {output.binny_dir} >> {log} 2>&1

        cd {output.binny_dir}
        echo "Changing config file." >> {log} 2>&1
 
        sed -i -e "s|prokka_env: \\"\\"|prokka_env: \\"{params.prokka_env}\\"|g" \
               -e "s|mantis_env: \\"\\"|mantis_env: \\"{params.mantis_env}\\"|g" \
               -e "s|conda_source: \\"\\"|conda_source: \\"{CONDA_DIR}\\"|g" \
               -e "s|snakemake_env: \\"\\"|snakemake_env: \\"{params.snakemake_env}\\"|g" config/config.*.yaml
        echo "Setting up binny." >> {log} 2>&1
        ./binny -i config/config.init.yaml >> {log} 2>&1 && touch {output.out}
        cd ../../
        """

##################################################
# SemiBin
##################################################

rule semibin_setup:
    input:
        "gtdbtk"
    output:
        "SemiBin_setup.done"
    log:
        os.path.join(INSTALL_DIR, "logs/SemiBin_setup.log")
    threads:
        2
    conda:
        os.path.join(MAIN_ENV_DIR, "IMP_binning.yaml")
    message:
        "SemiBin: Setting up SemiBin"
    shell:
        """
        echo TEST
        SemiBin \
            download_GTDB \
                --reference-db-data-dir {input} >> {log} 2>&1 && \
        SemiBin \
            check_install >> {log} 2>&1 && \
        touch {output}
        echo TEST2
        """

##################################################
# abawaca
##################################################

rule abawaca_setup:
    output:
        esom_dir=directory(os.path.join(SUBMOD_DIR, "ESOM")),
        out=os.path.join(SUBMOD_DIR, "ESOM/setup.done")
    log:
        os.path.join(INSTALL_DIR,"logs/abawaca_setup.log")
    threads:
        1
    conda:
        os.path.join(MAIN_ENV_DIR, "IMP_binning.yaml")
    message:
        "abawaca: Setting up ESOM for abawaca"
    shell:
        """
        echo "Downloading ESOM for abawaca." >> {log} 2>&1
        git clone https://github.com/tetramerFreqs/Binning {output.esom_dir} >> {log} 2>&1 && \
        touch {output.out}
        """
