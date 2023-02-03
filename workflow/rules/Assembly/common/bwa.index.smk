import math

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
        mem_mb = BIGMEMCORE_MB
    threads:
        # Ensure at least 20gb of mem
        max(int(math.ceil(20 / (int(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]]))), 1)
    params:
        fasta_file = lambda wc: "_".join(wc.get("fasta").split("/")[-2:]),
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "bwa_index: Indexing {params.fasta_file} for bwa."
    shell:
        """
        bwa index {wildcards.fasta}
        """
