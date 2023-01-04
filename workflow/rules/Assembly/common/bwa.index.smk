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
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/assembly_bwa_index_{fasta}.log"
    message: "bwa_index: Indexing {wildcards.fasta} for bwa."
    shell:
        """
        bwa index {wildcards.fasta} 
        """
