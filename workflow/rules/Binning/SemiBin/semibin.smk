rule semibin:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        mapping = "Assembly/mg.reads.sorted.bam"
    output:
        s2b = "Binning/SemiBin/contig2bin.tsv",
        out_dir = temp(directory("Binning/SemiBin")),
        out_bin_dir = directory("Binning/SemiBin/bins")
    params:
        model = config["binning"]["SemiBin"]["model"],
        cutoff = config["binning"]["SemiBin"]["cutoff"],
        # db_dir=os.path.join(config["db_path"], "gtdb_semibin"),
        tmp_dir = os.path.join(config["tmp_dir"])
    threads:
        getThreads(128)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        ENVDIR + "/IMP_semibin.yaml"
    log:
        "logs/binning/SemiBin.log"
    message:
        "Running SemiBin."
    shell:
        """
        # Run SemiBin2
        SemiBin2 \
            single_easy_bin \
                --input-fasta {input.assembly} \
                --input-bam {input.mapping} \
                --environment {params.model} \
                --threads {threads} \
                --random-seed 0 \
                --min-len {params.cutoff} \
                --compression none \
                --output {output.out_dir} >> {log} 2>&1

        # Move bin fastas and rename.
        mv {output.out_dir}/output_bins {output.out_dir}/bins
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
