rule maxbin_contig2bin:
    input:
        "Binning/MaxBin/maxbin_res.log"
    output:
        "Binning/MaxBin/maxbin_contig2bin.txt",
        "Binning/MaxBin/scaffold2bin.tsv"
    threads: 1
    params:
        runtime = "3:00:00",
        mem = MEMCORE
    log: "logs/binning_maxbin_contig2bin.log"
    message: "maxbin_contig2bin: Getting contig to bin mapping."
    shell:
        """
        for f in $(ls Binning/MaxBin/*.fasta | cut -f3 -d "/"); do 
         sed -e "s/$/\t$f/g" <(grep "^>" Binning/MaxBin/$f | sed -e 's/>//g') >> {output[0]}
        done
        ln -s ../../{output[0]} {output[1]}
        echo "Contig to bin mapping done" > {log}
        """
