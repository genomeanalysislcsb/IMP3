rule maxbin_contig2bin:
    input:
        "Binning/MaxBin/maxbin_res.log"
    output:
        "Binning/MaxBin/maxbin_contig2bin.txt",
        "Binning/MaxBin/scaffold2bin.tsv"
    threads: 1
    resources:
        runtime = "3:00:00",
        mem = MEMCORE
    log: "logs/binning_maxbin_contig2bin.log"
    message: "maxbin_contig2bin: Getting contig to bin mapping."
    shell:
        """
        for f in Binning/MaxBin/*.fasta; do
          if [ -s "$f" ]; then
            f2=`echo $f | cut -f3 -d "/"`
            sed -e "s/$/\t$f2/g" <(grep "^>" Binning/MaxBin/$f2 | sed -e 's/>//g') >> {output[0]}
          fi
        done
        if [ -e {output[0]} ]; then
            ln -s ../../{output[0]} {output[1]}
        else
            touch {output[0]} {output[1]}
        fi
        echo "Contig to bin mapping done" > {log}
        """
