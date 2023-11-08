rule metabat:
    input:
        assembly="Assembly/%s.assembly.merged.fa" % ASS,
        depth="Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
        s2b="Binning/MetaBAT/contig2bin.tsv",
        out_dir=directory("Binning/MetaBAT"),
        out_bin_dir=directory("Binning/MetaBAT/bins")
    params:
        tmp_dir=TMPDIR,
        res="Binning/MetaBAT/metabat_res",
        min_cont_size=config["binning"]["MetaBAT"]["cutoff"]
    threads:
        getThreads(8)
    resources:
        runtime = "04:00:00",
        mem_mb = MAIN_MEMCORE_MB
    conda:
        ENVDIR + "/IMP_metawrap.yaml"
    log:
        "logs/binning_MetaBAT.log"
    message:
        "Running MetaBAT."
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
