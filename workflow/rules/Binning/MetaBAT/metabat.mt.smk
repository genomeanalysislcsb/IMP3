rule metabat:
    input:
        assembly="Assembly/%s.assembly.merged.fa" % ASS
    output:
        s2b="Binning/MetaBAT/scaffold2bin.tsv",
        out_dir=directory("Binning/MetaBAT/bins")
    params:
        res="Binning/MetaBAT/metabat_res"
    threads:
        getThreads(8)
    resources:
        runtime = "04:00:00",
        mem_mb = MEMCORE_MB
    conda:
        ENVDIR + "/IMP_metawrap.yaml"
    log:
        "logs/binning_metabat.log"
    message:
        "metabat: Running MetaBAT."
    shell:
        """
        ## Create MetaBat dir
        TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "metabat_depth_XXXXXX.tsv")       
        ## Run MetaBat
        metabat2 -i {input.assembly} \
                 --saveCls \
                 -o {output.res} \
                 -t {threads} \
                 -m {config[binning][MetaBAT][cutoff]} > {log} 2>&1 
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
