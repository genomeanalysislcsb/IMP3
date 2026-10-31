
rule setup_binny:
    output:
        binny_ready_check = os.path.join(ROOTDIR, 'workflow/envs/binny.ready')
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        12
    conda:
        os.path.join(ENVDIR, "IMP_binny.yaml")
    log:
        os.path.join(ENVDIR, "logs/setup_binny.log")
    message:
        "binny: Setting up binny."
    shadow: "copy-minimal"
    shell:
        """
        # Setup binny
        echo "### Setting up binny ... ###" > {log}
        binny --setup >> {log} 2>&1

        # Touch ready file
        touch {output.binny_ready_check} >> {log} 2>&1

        echo "## Done. ###" >> {log}
        """

rule binny:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        depth = "Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
        s2b = "Binning/MetaBAT/contig2bin.tsv",
        out_dir = directory("Binning/MetaBAT"),
        out_bin_dir = directory("Binning/MetaBAT/bins")
    params:
        tmp_dir = TMPDIR,
        res = "Binning/MetaBAT/metabat_res",
        min_cont_size = config["binning"]["MetaBAT"]["cutoff"]
    threads:
        getThreads(8)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        ENVDIR + "/IMP_binny.yaml"
    log:
        "logs/binning/binny.log"
    message:
        "Running binny."
    shadow: "copy-minimal"
    shell:
        """
        ## Create MetaBat dir
        mkdir -p {output.out_bin_dir}

        export TMPDIR="{params.tmp_dir}"
        TMP_FILE=$(mktemp "metabat_depth_XXXXXX.tsv")
        echo "# fakeheader" > $TMP_FILE
        cat {input.depth}  >> $TMP_FILE

        ## Run MetaBat
        metabat2 -i {input.assembly} \
                 -a $TMP_FILE \
                 --cvExt \
                 --saveCls \
                 -o {params.res} \
                 -t {threads} \
                 -v \
                 -m {params.min_cont_size} >> {log} 2>&1
        mv {params.res} {output.s2b}
        rm $TMP_FILE
        
        # Move bin fastas and rename.
        shopt -s nullglob
        for i in {output.out_dir}/*.fa; do
          bin_name=$(basename -s .fa ${{i}})
          mv ${{i}} {output.out_bin_dir}/${{bin_name}}.fasta
        done
        shopt -u nullglob
        """

rule binny_contig2bin:
    input:
        binny_contig_data
    output:
        binny_s2b
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    message:
        "binny: Preparing table for DAS_Tool from binny output."
    run:
        import gzip
        with gzip.open(binny_contig_data,'rt') as gz: # {input.contig_data}
            with open(binny_s2b, 'w') as of:  # {output.s2b}
                header = next(gz)
                for line in gz:
                    line = line.strip('\n ')
                    line_contig = line.split('\t')[0]
                    line_bin = line.split('\t')[1]
                    if line_bin != 'N':
                        of.write(f'{line_contig}\t{line_bin}\n')