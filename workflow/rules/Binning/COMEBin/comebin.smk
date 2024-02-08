rule comebin:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        mapping = "Assembly/mg.reads.sorted.bam"
    output:
        s2b = "Binning/COMEBin/contig2bin.tsv",
        out_dir = temp(directory("Binning/COMEBin")),
        out_bin_dir = directory("Binning/COMEBin/bins")
    params:
        tmp_align_dir = "Binning/COMEBin/bam",
        abs_out_dir = OUTPUTDIR
    threads:
        getThreads(128)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        ENVDIR + "/IMP_comebin.yaml"
    log:
        os.path.join(OUTPUTDIR, "logs/binning/COMEBin.log")
    message:
        "Running COMEBin."
    shell:
        """
        # Change to outdir to avoid weird errors
        cd {params.abs_out_dir}/{output.out_dir}

        # Prepare tmp dir for mg read alignments
        mkdir -p {params.abs_out_dir}/{params.tmp_align_dir} >> {log} 2>&1
        ln -s {params.abs_out_dir}/{input.mapping} {params.abs_out_dir}/{params.tmp_align_dir}/ >> {log} 2>&1

        # Run COMEBin
        CUDA_VISIBLE_DEVICES=0 \
        run_comebin.sh -a {params.abs_out_dir}/{input.assembly} \
                       -o {params.abs_out_dir}/{output.out_dir} \
                       -p {params.abs_out_dir}/{params.tmp_align_dir} \
                       -t {threads} >> {log} 2>&1

        # Remove symlink dir
        rm -rf {params.abs_out_dir}/{params.tmp_align_dir} >> {log} 2>&1

        # Move bin fastas and rename.
        mv {params.abs_out_dir}/{output.out_dir}/comebin_res/comebin_res_bins {params.abs_out_dir}/{output.out_dir}/bins >> {log} 2>&1
        shopt -s nullglob
        for i in {params.abs_out_dir}/{output.out_bin_dir}/*.fa; do
          bin_name=$(basename -s .fa ${{i}})
          mv ${{i}} {params.abs_out_dir}/{output.out_bin_dir}/COMEBIN_${{bin_name}}.fasta >> {log} 2>&1
        done
        shopt -u nullglob

        # Write bin table
        python {SRCDIR}/write_bin_table_from_fastas.py {params.abs_out_dir}/{output.out_bin_dir} {params.abs_out_dir}/{output.s2b} >> {log} 2>&1
        rm -rf {params.abs_out_dir}/{output.out_dir}/comebin_res/comebin_res_bins >> {log} 2>&1
        """
