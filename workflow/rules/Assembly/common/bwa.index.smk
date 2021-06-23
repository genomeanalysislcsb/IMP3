rule bwa_index:
    input:
        "{fasta}"
    output:
        "{fasta}.amb",
        "{fasta}.bwt",
        "{fasta}.pac",
        "{fasta}.sa",
        "{fasta}.ann"
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/assembly_bwa_index.{fasta}.log"
    message: "bwa_index: Indexing {wildcards.fasta} for bwa."
    shell:
        """
        bwa index {wildcards.fasta} > {log} 2>&1
        """
