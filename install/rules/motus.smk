
rule motus_db_download:
    output:
        dl_done = "motus_db_download.done",
        db_dir = directory("motus_metadata"),
        fasta = "motus_metadata/db_mOTU_DB_CEN.fasta",
        coordinates = "motus_metadata/db_mOTU_padding_coordinates_CEN.tsv",
        metaphy = "motus_metadata/db_mOTU_taxonomy_meta-mOTUs.tsv",
        refphy = "motus_metadata/db_mOTU_taxonomy_ref-mOTUs.tsv"
    log:
        os.path.join(INSTALL_DIR, "logs/motus_db_download.log")
    threads:
        4
    params:
        lib_path = "lib/python3.9/site-packages/motus/db_mOTU",
        fasta = "db_mOTU_DB_CEN.fasta",
        coordinates = "db_mOTU_padding_coordinates_CEN.tsv",
        metaphy = "db_mOTU_taxonomy_meta-mOTUs.tsv",
        refphy = "db_mOTU_taxonomy_ref-mOTUs.tsv"
    conda:
        os.path.join(MAIN_ENV_DIR, "IMP_motus.yaml")
    message:
        "mOTUs: Downloading database."
    shell:
        """
        motus downloadDB > {log} 2>&1
        cp ${{CONDA_PREFIX}}/{params.lib_path}/{params.fasta} {output.fasta}
        cp ${{CONDA_PREFIX}}/{params.lib_path}/{params.coordinates} {output.coordinates}
        cp ${{CONDA_PREFIX}}/{params.lib_path}/{params.metaphy} {output.metaphy}
        cp ${{CONDA_PREFIX}}/{params.lib_path}/{params.refphy} {output.refphy}
        touch {output.dl_done}
        """
