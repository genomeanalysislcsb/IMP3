rule maxbin:
    input:
       assembly = "Assembly/%s.assembly.merged.fa" % ASS,
       depth = "Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
       "Binning/MaxBin/maxbin_res.log",
       "Binning/MaxBin/maxbin_res.summary",
       temp(directory("Binning/MaxBin"))
    threads:
        getThreads(14)
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        ENVDIR + "/IMP_metawrap.yaml"
    log:
        "logs/binning/maxbin.log"
    message:
        "maxbin: Running MaxBin2."
    shell:
        """
        export PATH=${{CONDA_PREFIX}}/bin:$PATH
        export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl
        # Run MaxBin
        run_MaxBin.pl -contig {input.assembly} \
                      -abund {input.depth} \
                      -out Binning/MaxBin/maxbin_res \
                      -thread {threads} \
                      -min_contig_length {config[binning][MaxBin][cutoff]} > {log} 2>&1 ||\
                       touch {output} status/binning_maxbin.impossible
        """
