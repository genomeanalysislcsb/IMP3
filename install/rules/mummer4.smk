rule install_qm_mummer4:
    output:
        "qm.ready"
    log:
        os.path.join(INSTALL_DIR, "logs", "quickmerge_install.log")
    threads:
        1
    conda:
        os.path.join(MAIN_ENV_DIR, "IMP_assembly_merge.yaml")
    message:
        "Fixing quickmerge wrapper."
    shell:
        """
        wget https://raw.githubusercontent.com/mahulchak/quickmerge/master/merge_wrapper.py \
             -O ${{CONDA_PREFIX}}/bin/merge_wrapper.py >> {log} 2>&1 \
          && chmod 750 ${{CONDA_PREFIX}}/bin/merge_wrapper.py \
          && touch {output}
        """