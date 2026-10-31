rule maxbin_contig2bin:
    input:
        "Binning/MaxBin/maxbin_res.log"
    output:
        c2b="Binning/MaxBin/maxbin_contig2bin.txt",
        s2b="Binning/MaxBin/contig2bin.tsv",
        out_dir=directory('Binning/MaxBin/bins')
    threads: 1
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log: "logs/binning/maxbin_contig2bin.log"
    message: "maxbin_contig2bin: Getting contig to bin mapping."
    shadow: "copy-minimal"
    shell:
        """
        for f in Binning/MaxBin/*.fasta; do
          if [ -s "$f" ]; then
            f2=`echo $f | cut -f3 -d "/"`
            sed -e "s/$/\t$f2/g" <(grep "^>" Binning/MaxBin/$f2 | sed -e 's/>//g') >> {output.c2b}
          fi
        done
        if [ -e {output[0]} ]; then
            ln -s ../../{output.c2b} {output.s2b}
        else
            touch {output.c2b} {output.s2b}
        fi
        echo "Contig to bin mapping done" > {log}
        mkdir -p {output.out_dir} && mv Binning/MaxBin/*.fasta {output.out_dir}
        """
