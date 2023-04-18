rule tar_intermediary:
    input:
        # Make sure its done at the very end
        "status/workflow.done",
        "Assembly/%s.assembly.merged.fa" % ASS
    output:
        "Assembly/intermediary.tar.gz"
    threads:
        1
    resources:
        runtime = "04:00:00",
        mem_mb = MEMCORE_MB
    params:
        intermediary = "Assembly/intermediary/"
    log:
        "logs/assembly_tar_intermediary.log"
    message:
        "tar_intermediary: Compressing intermediary steps of assembly."
    shell:
       """
       tar czf {output} {params.intermediary} && rm -r {params.intermediary} >> {log} 2>&1
       """

rule contig_fasta_indexing:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS
    output:
        "Assembly/%s.assembly.merged.fa.fai" % ASS
    resources:
        runtime = "04:00:00",
        mem_mb = MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/assembly_contig_fasta_indexing.log"
    message:
        "contig_fasta_indexing: Indexing assembly."
    shell:
        """
        samtools faidx {input[0]} > {log} 2>&1
        """

rule contig_fasta2bed_conversion:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS,
        "Assembly/%s.assembly.merged.fa.fai" % ASS
    output:
        "Assembly/%s.assembly.merged.fa.bed3" % ASS
    resources:
        runtime = "02:00:00",
        mem_mb = MEMCORE_MB
    threads:
        1
    # log:
    #     "logs/assembly_fastabed.log"
    message:
        "contig_fasta2bed_conversion: Writing bed file for assembly."
    shell:
        """
        cat {input[0]}.fai | awk '{{print $1 \"\t0\t\" $2}}' > {output}
        """
