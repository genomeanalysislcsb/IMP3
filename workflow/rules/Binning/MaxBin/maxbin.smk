rule maxbin:
    input:
       assembly = "Assembly/%s.assembly.merged.fa" % ASS,
       depth = "Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
        mb_summary = "Binning/MaxBin/maxbin_res.summary",
        out_dir = temp(directory("Binning/MaxBin")),
        c2b = "Binning/MaxBin/maxbin_contig2bin.txt",
        s2b = "Binning/MaxBin/contig2bin.tsv",
        out_bin_dir = directory('Binning/MaxBin/bins')
    threads:
        getThreads(14)
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        min_cont_len = config['binning']['MaxBin']['cutoff']
    conda:
        ENVDIR + "/IMP_maxbin.yaml"
    log:
        "logs/binning/maxbin.log"
    message:
        "maxbin: Running MaxBin2."
    shell:
        """
        # Export env vars
        export PATH=${{CONDA_PREFIX}}/bin:$PATH 2> {log}
        export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl 2> {log}

        # Run MaxBin
        echo "Running MaxBin2 ..." >> {log}
        run_MaxBin.pl -contig {input.assembly} \
                      -abund {input.depth} \
                      -out {output.out_dir}/maxbin_res \
                      -thread {threads} \
                      -min_contig_length {params.min_cont_len} >> {log} 2>&1 \
          || (touch {output.mb_summary} \
                    {output.c2b} \
                    {output.s2b} \
                    status/binning_maxbin.impossible >> {log} 2>&1; \
              mkdir -p {output.out_dir} \
                       {output.out_bin_dir} >> {log} 2>&1
             ) \
          && (echo "Cannot run MaxBin2" >> {log}" && exit 0)

        for f in {output.out_dir}/*.fasta; do
          if [ -s "$f" ]; then
            f2=`echo $f | cut -f3 -d "/"`
            sed -e "s/$/\t$f2/g" <(grep "^>" {output.out_dir}/$f2 | sed -e 's/>//g') >> {output.c2b}
          fi
        done

        echo "Contig to bin mapping done" >> {log}
        mkdir -p {output.out_bin_dir} && mv {output.out_dir}/*.fasta {output.out_bin_dir}
        
        echo "MaxBin2 done" >> {log}
        """
