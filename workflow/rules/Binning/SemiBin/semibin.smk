rule semibin:
    input:
        assembly="Assembly/%s.assembly.merged.fa" % ASS,
        mapping="Assembly/mg.reads.sorted.bam"
    output:
        s2b="Binning/SemiBin/scaffold2bin.tsv",
        out_dir=directory("Binning/SemiBin"),
        out_bin_dir=directory("Binning/SemiBin/bins")
    params:
#         model=config["binning"]["SemiBin"]["model"],
        cutoff=config["binning"]["SemiBin"]["cutoff"],
        db_dir=os.path.join(config["db_path"], "gtdb_semibin"),
        tmp_dir=os.path.join(config["tmp_dir"])
    threads:
        getThreads(128)
    resources:
        runtime = "12:00:00",
        mem_mb = MEMCORE_MB
    conda:
        ENVDIR + "/IMP_semibin.yaml"
    log:
        "logs/binning_SemiBin.log"
    message:
        "Running SemiBin."
    shell:
        """
        # Run SemiBin
        SemiBin \
            single_easy_bin \
                --input-fasta {input.assembly} \
                --input-bam {input.mapping} \
                --self-supervised \
                --threads {threads} \
                --random-seed 0 \
                --min-len {params.cutoff} \
                --reference-db-data-dir {params.db_dir} \
                --output {output.out_dir} >> {log} 2>&1

        # Move bin fastas and rename.
        mv {output.out_dir}/output_recluster_bins {output.out_dir}/bins
        shopt -s nullglob
        for i in {output.out_bin_dir}/*.fa; do
          bin_name=$(basename -s .fa ${{i}})
          mv ${{i}} {output.out_bin_dir}/${{bin_name}}.fasta
        done
        shopt -u nullglob
        # Write bin table
        python {SRCDIR}/write_bin_table_from_fastas.py {output.out_bin_dir} {output.s2b}
        rm -rf {output.out_dir}/output_bins
        """
