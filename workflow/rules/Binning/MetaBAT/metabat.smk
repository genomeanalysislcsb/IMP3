rule metabat:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        depth = "Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
        s2b = "Binning/MetaBAT/contig2bin.tsv",
        out_dir = directory("Binning/MetaBAT"),
        out_bin_dir = directory("Binning/MetaBAT/bins"),
        ass_stats = temp("Binning/MetaBAT/assembly_stats.tsv")
    params:
        tmp_dir = TMPDIR,
        res = "Binning/MetaBAT/metabat_res",
        min_cont_size = config["binning"]["MetaBAT"]["cutoff"],
        impossible_check = f'{OUTPUTDIR}/status/binning_metabat.impossible'
    threads:
        getThreads(8)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        ENVDIR + "/IMP_metabat.yaml"
    log:
        "logs/binning/MetaBAT.log"
    message:
        "Running MetaBAT."
    shadow: "copy-minimal"
    shell:
        """
        # Clear previous status check
        rm -rf {params.impossible_check} > {log} 2>&1

        # Check if there are enough contigs < min_cont_size to proceed with binning
        echo "### Checking number of contigs > min length ... ###" >> {log}
        # Stats format is `file  format  type  num_seqs    sum_len  min_len  avg_len  max_len`
        seqkit seq -m {params.min_cont_size} {input.assembly} 2>> {log} | seqkit stats -j 4 -T 2>> {log} | tail -n +2 > {output.ass_stats} 2>> {log}
        ass_n_contigs=$(cut -f 4 {output.ass_stats} 2>> {log})
        ass_size=$(cut -f 5 {output.ass_stats} 2>> {log})

        # Touch output and extit if theres less than 1000 contigs or the size is < 2*10^6 bps
        if [[ ! $ass_n_contigs -ge 1000 || ! $ass_size -ge 2000000 ]];then
          echo "### Less than 1000 contigs larger min size or assembly smaller than 2*10^6 bps. Exiting ... ###" >> {log}
          touch {output.s2b} \
                {params.impossible_check} >> {log} 2>&1
          mkdir -p {output.out_dir} \
                   {output.out_bin_dir} >> {log} 2>&1
          exit 0
        fi

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
