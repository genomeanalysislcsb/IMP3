rule metabat:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS,
        "Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
        "Binning/MetaBAT/metabat_res",
	"Binning/MetaBAT/scaffold2bin.tsv"
    threads: 4
    params:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    conda: ENVDIR + "/IMP_binning.yaml"
    log: "logs/binning_MetaBAT.log"
    message: "Running MetaBAT."
    shell:
        """
        ## Create MetaBat dir
        mkdir -p Binning/MetaBAT
	TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "metabat_depth_XXXXXX.tsv")       
 	echo "# fakeheader" > $TMP_FILE
	cat {input[1]}  >> $TMP_FILE
        ## Run MetaBat
        metabat2 -i {input[0]} \
         -a $TMP_FILE --cvExt \
	 --saveCls \
         -o {output[0]} \
         -t {threads} -v \
         -m {config[binning][MetaBAT][cutoff]} > {log} 2>&1
	ln -s ../../{output[0]} {output[1]}
	rm $TMP_FILE
        """
