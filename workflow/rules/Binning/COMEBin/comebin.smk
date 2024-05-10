# TODO: Manually changed comebin cluster.py to try to deal with cases where no seed clusters are found

# rule comebin:
#     input:
#         assembly = "Assembly/%s.assembly.merged.fa" % ASS,
#         mapping = "Assembly/mg.reads.sorted.bam"
#     output:
#         s2b = "Binning/COMEBin/contig2bin.tsv",
#         out_dir = directory("Binning/COMEBin"),
#         out_bin_dir = directory("Binning/COMEBin/bins"),
#         filt_ass = temp(f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa')
#     params:
#         tmp_align_dir = "Binning/COMEBin/bam",
#         abs_out_dir = OUTPUTDIR,
#         min_cont_size = config["binning"]["comebin"]["cutoff"],
#         # filt_ass = f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa'
#     threads:
#         # getThreads(7, high_mem_cluster=True)
#         getThreads(128)
#     resources:
#         runtime = LONG_RUNTIME,
#         # mem_mb = 27000,
#         mem_mb = MAIN_MEMCORE_MB
#     conda:
#         ENVDIR + "/IMP_comebin.yaml"
#     log:
#         os.path.join(OUTPUTDIR, "logs/binning/COMEBin.log")
#     message:
#         "Running COMEBin."
#     shell:
#         """
#         # Change to outdir to avoid weird errors
#         echo "Moving to {params.abs_out_dir}/{output.out_dir} ..." >> {log} 2>&1 
#         cd {params.abs_out_dir}/{output.out_dir} >> {log} 2>&1

#         # Prepare tmp dir for mg read alignments
#         echo "Symlinking mg read alignments to tmp dir {params.abs_out_dir}/{params.tmp_align_dir} ..." >> {log} 2>&1
#         mkdir -p {params.abs_out_dir}/{params.tmp_align_dir} >> {log} 2>&1
#         ln -s {params.abs_out_dir}/{input.mapping} {params.abs_out_dir}/{params.tmp_align_dir}/ >> {log} 2>&1

#         # Filter contigs by length
#         echo "Filtering contigs by length ..." >> {log} 2>&1
#         $CONDA_PREFIX/bin/COMEBin/scripts/Filter_tooshort.py {params.abs_out_dir}/{input.assembly} {params.min_cont_size} >> {log} 2>&1
#         echo "Filtered assembly written to {params.abs_out_dir}/{output.filt_ass}." >> {log} 2>&1

#         # Run COMEBin
#         echo "Running COMEBin ..." >> {log} 2>&1
#         CUDA_VISIBLE_DEVICES=0

#         run_comebin.sh -a {params.abs_out_dir}/{output.filt_ass} \
#                        -o {params.abs_out_dir}/{output.out_dir} \
#                        -p {params.abs_out_dir}/{params.tmp_align_dir} \
#                        -t {threads} >> {log} 2>&1

#         # Remove symlink dir
#         echo "Removing tmp read map dir ..." >> {log} 2>&1
#         rm -rf {params.abs_out_dir}/{params.tmp_align_dir} >> {log} 2>&1

#         # Move bin fastas and rename.
#         echo "Moving and renaming bins ..." >> {log} 2>&1
#         mv {params.abs_out_dir}/{output.out_dir}/comebin_res/comebin_res_bins {params.abs_out_dir}/{output.out_dir}/bins >> {log} 2>&1
#         shopt -s nullglob
#         for i in {params.abs_out_dir}/{output.out_bin_dir}/*.fa; do
#           bin_name=$(basename -s .fa ${{i}})
#           mv ${{i}} {params.abs_out_dir}/{output.out_bin_dir}/COMEBIN_${{bin_name}}.fasta >> {log} 2>&1
#         done
#         shopt -u nullglob

#         # Write bin table
#         echo "Writing bin table ..." >> {log} 2>&1
#         python {SRCDIR}/write_bin_table_from_fastas.py {params.abs_out_dir}/{output.out_bin_dir} {params.abs_out_dir}/{output.s2b} >> {log} 2>&1
#         rm -rf {params.abs_out_dir}/{output.out_dir}/comebin_res/comebin_res_bins >> {log} 2>&1
#         echo "Done." >> {log} 2>&1
#         """

rule comebin_train:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        mapping = "Assembly/mg.reads.sorted.bam"
    output:
        out_dir = directory("Binning/COMEBin"),
        filt_ass = temp(f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa'),
        filt_ass_lst = temp(f'Assembly/contig_length_filter1000.txt'),
        filt_ass_len = temp(f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa_lengths.txt'),
        tmp_align_dir = temp(directory('Binning/COMEBin/bam')),
        tmp_feature_dir = temp(directory('Binning/COMEBin/data_augmentation')),
        tmp_res_dir = temp(directory('Binning/COMEBin/comebin_res')),
        embeddings = temp('Binning/COMEBin/comebin_res/embeddings.tsv'),
    params:
        abs_out_dir = OUTPUTDIR,
        min_cont_size = config["binning"]["comebin"]["cutoff"]
    threads:
        getThreads(7, high_mem_cluster=True)
    resources:
        runtime = 480,
        mem_mb = 27000
    conda:
        os.path.join(ENVDIR, "IMP_comebin.yaml")
    log:
        os.path.join(OUTPUTDIR, "logs/binning/COMEBin_train.log")
    message:
        "Running COMEBin training."
    shell:
        """
        # Change to outdir to avoid weird errors
        echo "Moving to {params.abs_out_dir}/{output.out_dir} ..." >> {log} 2>&1 
        cd {params.abs_out_dir}/{output.out_dir} >> {log} 2>&1

        # Prepare tmp dir for mg read alignments
        echo "Symlinking mg read alignments to tmp dir {params.abs_out_dir}/{output.tmp_align_dir} ..." >> {log} 2>&1
        mkdir -p {params.abs_out_dir}/{output.tmp_align_dir} >> {log} 2>&1
        ln -s {params.abs_out_dir}/{input.mapping} {params.abs_out_dir}/{output.tmp_align_dir}/ >> {log} 2>&1

        # Filter contigs by length
        echo "Filtering contigs by length ..." >> {log} 2>&1
        $CONDA_PREFIX/bin/COMEBin/scripts/Filter_tooshort.py {params.abs_out_dir}/{input.assembly} {params.min_cont_size} >> {log} 2>&1
        echo "Filtered assembly written to {params.abs_out_dir}/{output.filt_ass}." >> {log} 2>&1

        # Run COMEBin
        echo "Running COMEBin training ..." >> {log} 2>&1

        run_comebin_train.sh -a {params.abs_out_dir}/{output.filt_ass} \
                             -o {params.abs_out_dir}/{output.out_dir} \
                             -p {params.abs_out_dir}/{output.tmp_align_dir} \
                             -t {threads} >> {log} 2>&1

        echo "Done." >> {log} 2>&1
        """

rule comebin_bin:
    input:
        out_dir = 'Binning/COMEBin',
        filt_ass = f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa',
        filt_ass_lst = f'Assembly/contig_length_filter1000.txt',
        filt_ass_len = f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa_lengths.txt',
        tmp_feature_dir = 'Binning/COMEBin/data_augmentation',
        tmp_res_dir = 'Binning/COMEBin/comebin_res',
        embeddings = 'Binning/COMEBin/comebin_res/embeddings.tsv',
    output:
        s2b = "Binning/COMEBin/contig2bin.tsv",
        out_bin_dir = directory("Binning/COMEBin/bins"),
        seed_file = temp(f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa.bacar_marker.2quarter_lencutoff_1001.seed'),
    params:
        abs_out_dir = OUTPUTDIR
    threads:
        getThreads(128)
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_comebin.yaml")
    log:
        os.path.join(OUTPUTDIR, "logs/binning/COMEBin_bin.log")
    message:
        "Running COMEBin binning."
    shell:
        """
        # Change to outdir to avoid weird errors
        echo "Moving to {params.abs_out_dir}/{input.out_dir} ..." >> {log} 2>&1 
        cd {params.abs_out_dir}/{input.out_dir} >> {log} 2>&1

        # Run COMEBin binning
        echo "Running COMEBin binning ..." >> {log} 2>&1
        CUDA_VISIBLE_DEVICES=0

        run_comebin_bin.sh {params.abs_out_dir}/{input.out_dir} \
                           {params.abs_out_dir}/{input.filt_ass} \
                           {threads} >> {log} 2>&1

        # Move bin fastas and rename.
        echo "Moving and renaming bins ..." >> {log} 2>&1
        mv {params.abs_out_dir}/{input.out_dir}/comebin_res/comebin_res_bins {params.abs_out_dir}/{input.out_dir}/bins >> {log} 2>&1
        shopt -s nullglob
        for i in {params.abs_out_dir}/{output.out_bin_dir}/*.fa; do
          bin_name=$(basename -s .fa ${{i}})
          mv ${{i}} {params.abs_out_dir}/{output.out_bin_dir}/COMEBIN_${{bin_name}}.fasta >> {log} 2>&1
        done
        shopt -u nullglob

        # Write bin table
        echo "Writing bin table ..." >> {log} 2>&1
        python {SRCDIR}/write_bin_table_from_fastas.py {params.abs_out_dir}/{output.out_bin_dir} {params.abs_out_dir}/{output.s2b} >> {log} 2>&1
        
        echo "Done." >> {log} 2>&1
        """
