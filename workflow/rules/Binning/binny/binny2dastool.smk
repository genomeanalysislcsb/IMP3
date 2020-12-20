
rule binny_pgol_table:
    input:
        expand("Binning/binny/contigs2clusters.{pk}.{nn}.tsv",pk=config["binning"]["binny"]["pk"],nn=config["binning"]["binny"]["nn"]),
        "Binning/binny/clusterFiles.tar.gz"
    output:
        "Binning/binny/scaffold2bin.tsv"
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    threads: 1
#    log: "logs/binning_binny.log"
    message: "binny_pgol_table: Preparing table for DAStool from binny output."
    shell:
        """
        grep [PGOL] {input[0]} >> {output} || touch {output}
        """

