# rule get_ccp_crap_db:
#     output:
#         temp("Analysis/annotation/latest_CCP_cRAP.fasta")
#     conda:
#         os.path.join(ENVDIR, "IMP_proteomics.yaml")
#     log:
#         "logs/analysis/get_ccp_crap_db.log"
#     message:
#         "Getting the CCP cRAP database in R."
#     shadow: "copy-minimal"
#     shell:
#         """
#         Rscript -e 'remotes::install_github("CambridgeCentreForProteomics/camprotR", dependencies=TRUE)' \
#                 -e 'library(Biostrings); library(camprotR); library(httr); download_ccp_crap("2023_01_CCP_cRAP.fasta", is_crap=TRUE, verbose=TRUE)' \
#         > {log} 2>&1
#         """

rule get_sipros_ensemble:
    output:
        temp("Analysis/annotation/Sipros-Ensemble")
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    log:
        "logs/analysis/get_sipros_ensemble.log"
    message:
        "Getting the Sipros-Ensemble tool."
    shadow: "copy-minimal"
    shell:
        """
        git clone https://github.com/ohickl/Sipros-Ensemble.git
        """

rule sipros_prepare_db:
    input:
        db_raw = "<some_path>/all_sample_db/mibipa_all_hs_crap_final_num_header.faa",
        config = "<some_path>/1229_2_sipros_test_all_samp_db_01/config/1229_2_sipros.cfg"
    output:
        db = "<some_path>/1229_2_sipros_test_all_samp_db_01/db/mibipa_all_hs_crap_final_num_header_sipros_rev.faa"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        py2 = os.path.join(CONDA_PREFIX, "bin/python2.7"),
        sip_src = "/mnt/data/local_tools/Sipros-Ensemble/Scripts"
    threads:
        1
    conda:
        os.path.join(ENVDIR, "python2_environment.yaml")
    log:
        "logs/analysis/sipros_prepare_db.log"
    message:
        "Preparing the sipros database."
    shadow: "copy-minimal"
    shell:
        """
        {params.py2} {params.sip_src}/sipros_prepare_protein_database.py -i {input.db_raw} \
                                                                         -o {output.db} \
                                                                         -c {input.config} \
        > {log} 2>&1
        """

rule sipros:
    input:
        config = "<some_path>/1229_2_sipros_test_all_samp_db_01/config/1229_2_sipros.cfg",
        ms_data = directory("<some_path>/1229_2_sipros_test_01/data"),
        db = rules.sipros_prepare_db.output.db
    output:
        out_dir = directory("<some_path>/1229_2_sipros_test_all_samp_db_01/out")
    resources:
        runtime = LONG_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "sipros.yaml")
    log:
        "logs/analysis/run_sipros.log"
    message:
        "Running sipros."
    shadow: "copy-minimal"
    shell:
        """
        Sipros_OpenMP -w {input.ms_data} \
                      -c {input.config} \
                      -o {output.out_dir} \
        > {log} 2>&1
        """

rule run_sipros_filtering:
    input:
        config = "<some_path>/1229_2_sipros_test_all_samp_db_01/config/1229_2_sipros.cfg",
        out_dir = rules.run_sipros.output.out_dir
    output:
        filtered_out = directory("<some_path>/1229_2_sipros_test_all_samp_db_01/out/filtered")
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "sipros_environment.yaml")
    log:
        "logs/analysis/run_sipros_filtering.log"
    message:
        "Running sipros filtering."
    shadow: "copy-minimal"
    shell:
        """
        runSiprosFiltering.sh -in {input.out_dir} \
                              -o {output.filtered_out} \
                              -c {input.config} \
        > {log} 2>&1
        """
