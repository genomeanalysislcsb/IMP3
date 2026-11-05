##################################################
# Checkm setup for MetaWrap
##################################################

rule setup_checkm:
    output:
        checkm_db=directory("checkm"),
        setup_check="checkm/checkm_db.ready"
    params:
        str1_to_fix='best_bin_set=binsA',
        str1_replacement='mv binsA binsO; mv binsA.stats binsO.stats; best_bin_set=binsO; rm -rf tmp.fa',
        str2_to_fix='\[ \$dereplicate \!\= \"false\" \]\; then',
        str2_replacement='\[ \$dereplicate \!\= \"false\" \] \&\& \[\[ \! \$n_binnings \-eq 1 \]\]; then',
        str3_to_fix='if \[\[ \$n_binnings \-ne 1 \]\]\; then',
        str3_replacement='if [[ \$n_binnings -ge 1 ]]; then',
        str4_to_fix='\${{bins1\#\#\*\/}}',
        str4_replacement='binsA',
        str5_to_fix='\${{bins2\#\#\*\/}}',
        str5_replacement='binsB',
        str6_to_fix='\${{bins3\#\#\*\/}}',
        str6_replacement='binsC'
    conda:
        os.path.join(MAIN_ENV_DIR, "IMP_metawrap.yaml")
    log:
        os.path.join(INSTALL_DIR, "logs/Metawrap_CheckM_setup.log")
    message:
        "MetaWrap: Setting up CheckM."
    shadow: "copy-minimal"
    shell:
        """
        echo "Downloading CheckM database" >> {log} 2>&1
        wget -nv https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -P {output.checkm_db} >> {log} 2>&1 && \
        echo "Decompressing CheckM database"  >> {log} 2>&1 && \
        tar -xvf {output.checkm_db}/checkm_data_2015_01_16.tar.gz -C {output.checkm_db} >> {log} 2>&1 && \
        rm -rf {output.checkm_db}/checkm_data_2015_01_16.tar.gz && \

        echo "Setting up CheckM" >> {log} 2>&1
        unset PYTHONHOME
        SETUPTOOLS_USE_DISTUTILS=stdlib
        printf "%s\n" {output.checkm_db} | checkm data setRoot {output.checkm_db} >> {log} 2>&1 && touch {output.setup_check}
        echo "Setup done." >> {log} 2>&1
        
        # Fix lines that lead to crash with MW, if theres only one binner output
        sed -i -e "331s/{params.str1_to_fix}/{params.str1_replacement}/g" \
               -e "397s/{params.str2_to_fix}/{params.str2_replacement}/g" \
               -e "451s/{params.str3_to_fix}/{params.str3_replacement}/g" \
               -e "462,463s/{params.str4_to_fix}/{params.str4_replacement}/g" \
               -e "465,466s/{params.str5_to_fix}/{params.str5_replacement}/g" \
               -e "469,470s/{params.str6_to_fix}/{params.str6_replacement}/g" \
               ${{CONDA_PREFIX}}/bin/metawrap-modules/bin_refinement.sh
        
        # Adjust score formula from 'completeness - (contamination * 5)' to 'completeness - (contamination * 2)'
        sed -i -e "s/\*5/*2/g" ${{CONDA_PREFIX}}/bin/metawrap-scripts/consolidate_two_sets_of_bins.py
        """

##################################################
# binny
##################################################
import errno

sm_user_path = config["binny"]["snakemake_env_path"]

if os.path.exists(sm_user_path) or sm_user_path == "in_path":
    sm_env_path = sm_user_path
else:
    raise FileNotFoundError(
    errno.ENOENT, os.strerror(errno.ENOENT), sm_user_path)

rule binny_setup:
    input:
        "mantis/setup.done"
    output:
        binny_dir=directory(os.path.join(SUBMOD_DIR, "binny")),
        out=os.path.join(SUBMOD_DIR, "binny/setup.done")
    log:
        os.path.join(INSTALL_DIR, "logs/binny_setup.log")
    params:
        prokka_env=os.path.join(MAIN_ENV_DIR,"IMP_annotation.yaml"),
        mantis_env=os.path.join(MAIN_ENV_DIR, "IMP_mantis.yaml"),
        snakemake_env=sm_env_path
    threads:
        4
    message:
        "binny: Setting up binny"
    shadow: "copy-minimal"
    shell:
        """
        echo "Downloading binny." >> {log} 2>&1
        git clone --quiet -b imp https://github.com/ohickl/binny.git \
        {output.binny_dir} >> {log} 2>&1

        cd {output.binny_dir}

        # Prevent binny subworkflow from changing IMP out dir
        sed -i -e "s|workdir\:|# workdir:|g" \
               -e "s|    OUTPUTDIR$|#     OUTPUTDIR|g" Snakefile
        
        echo "Changing config file." >> {log} 2>&1
 
        sed -i -e "s|^prokka_env: \\"\\"|prokka_env: \\"{params.prokka_env}\\"|g" \
               -e "s|^mantis_env: \\"\\"|mantis_env: \\"{params.mantis_env}\\"|g" \
               -e "s|^conda_source: \\"\\"|conda_source: \\"{CONDA_DIR}\\"|g" \
               -e "s|^snakemake_env: \\"\\"|snakemake_env: \\"{params.snakemake_env}\\"|g" config/config.*.yaml
        echo "Setting up binny." >> {log} 2>&1
        ./binny -i config/config.init.yaml >> {log} 2>&1 && touch {output.out}
        cd ../../
        """

##################################################
# SemiBin
##################################################

rule semibin_setup:
    output:
        gtdb_dir=directory("gtdb_semibin"),
        semibin_install_done="SemiBin_setup.done"
    log:
        os.path.join(INSTALL_DIR, "logs/SemiBin_setup.log")
    threads:
        2
    conda:
        os.path.join(MAIN_ENV_DIR, "IMP_semibin.yaml")
    message:
        "SemiBin: Setting up SemiBin"
    shadow: "copy-minimal"
    shell:
        """
        mkdir -p {output.gtdb_dir}
        SemiBin \
            download_GTDB \
                --reference-db-data-dir {output.gtdb_dir} >> {log} 2>&1 && \
        SemiBin \
            check_install >> {log} 2>&1 && \
        touch {output.semibin_install_done}
        """
