rule bwa_index:
    input:
        "{fasta}"
    output:
        "{fasta}.bwt.2bit.64",
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_mapping.yaml"
#    log: "logs/assembly_bwa_index.{fasta}.log"
    message: "bwa_index: Indexing {wildcards.fasta} for bwa-mem2."
    shell:
        """
        bwa-mem2 index {wildcards.fasta}
        """
