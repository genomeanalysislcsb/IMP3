import math

rule bwa_index:
    input:
        "{fasta}"
    output:
        "{fasta}.0123",
        "{fasta}.amb",
        "{fasta}.ann",
        "{fasta}.bwt.2bit.64",
        "{fasta}.pac"
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        # Ensure at least 100gb of mem
        max(int(math.ceil(100 / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    params:
        fasta_file = lambda wildcards: "_".join(wildcards.get("fasta").split("/")[-2:]),
        workdir = OUTPUTDIR
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "bwa_index: Indexing {params.fasta_file} for bwa."
    shell:
        """
        log="{params.workdir}/logs/preprocessing/index_{params.fasta_file}.log"

        # Ensure log dir exists
        mkdir -p {params.workdir}/logs/preprocessing

        echo "[$(date)] Start indexing {params.fasta_file} for bwa." >> $log
        bwa-mem2 index {wildcards.fasta} >> $log 2>&1
        echo "[$(date)] Finished indexing {params.fasta_file} for bwa." >> $log
        """
