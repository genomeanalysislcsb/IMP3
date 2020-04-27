rule maxbin:
    input:
       "Assembly/%s.assembly.merged.fa" % ASS,
       "Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
       "Binning/MaxBin/maxbin_res.log",
       "Binning/MaxBin/maxbin_res.summary"
    threads: 8
    params:
        runtime = "24:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_binning.yaml"
    log: "logs/binning_maxbin.log"
    message: "maxbin: Running MaxBin."
    shell:
        """
	## Create MaxBin dir
#        mkdir -p Binning/MaxBin
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
	## Run MaxBin
        run_MaxBin.pl -contig {input[0]} \
         -abund {input[1]} \
         -out Binning/MaxBin/maxbin_res \
         -thread {threads} \
         -min_contig_length {config[binning][MaxBin][cutoff]} > {log} 2>&1 || touch {output} status/binning_maxbin.impossible
        """
