rule tar_assembly_intermediary:
    input:
        # Make sure its done at the very end
        check = "status/workflow.done",
        assembly = "Assembly/%s.assembly.merged.fa" % ASS
    output:
        "Assembly/intermediary.tar.gz"
    threads:
        getThreads(128, bigmem=False, high_mem_cluster=False)
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        intermediary = "Assembly/intermediary/"
    log:
        "logs/assembly/tar_assembly_intermediary.log"
    message:
        "tar_intermediary: Compressing intermediary steps of assembly."
    shell:
        """
        echo "Compressing assembly intermediary dir ..." >> {log}

        tar -cf - -C Assembly/ intermediary/  2>> {log} | pigz -p {threads} -c > {output} 2>> {log} \
            && rm -rf {params.intermediary} >> {log} 2>&1

        echo "Done." >> {log}
        """

rule contig_fasta_indexing:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS
    output:
        "Assembly/%s.assembly.merged.fa.fai" % ASS
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/assembly/contig_fasta_indexing.log"
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
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    log:
        "logs/assembly/convert_fasta2bed.log"
    message:
        "contig_fasta2bed_conversion: Writing bed file for assembly."
    shell:
        """
        cat {input[0]}.fai | awk '{{print $1 \"\t0\t\" $2}}' > {output} 2> {log}
        """
