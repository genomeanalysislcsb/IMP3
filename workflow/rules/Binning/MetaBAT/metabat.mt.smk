rule metabat:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS
    output:
        "Binning/MetaBAT/metabat_res",
	"Binning/MetaBAT/scaffold2bin.tsv"
    threads: getThreads(8)
    resources:
        runtime = "24:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_binning.yaml"
    log: "logs/binning_metabat.log"
    message: "metabat: Running MetaBAT."
    shell:
        """
        ## Create MetaBat dir
	TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "metabat_depth_XXXXXX.tsv")       
        ## Run MetaBat
        metabat2 -i {input} \
	 --saveCls \
         -o {output[0]} \
         -t {threads} \
         -m {config[binning][MetaBAT][cutoff]} > {log} 2>&1 
	ln -s ../../{output[0]} {output[1]}
	rm $TMP_FILE
        """
