rule preproc_filter_genomes:
    input:
        reads=filtering_reads_se_input,
        ref=filtering_filter
    output:
        "Preprocessing/{omic}.se.{filterstep}.fq"
    log:
        "logs/preproc_filter_genomes.{omic}.{filterstep}.log"
    benchmark:
        "benchmark/preproc_filter_genomes.{omic}.{filterstep}.tsv"
    wildcard_constraints:
        omic="|".join(get_omics()),
        filterstep="|".join(["trimmed(\.rna_filtered)?\." + s for s in get_read_files_filtered_genomes_basenames(True)]),
    threads: get_cores_big()
    resources:
        runtime="08:00:00",
        mem=get_mem_big,
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "Filter reads: {input.reads} against {input.ref}"
    shell:
        """
        bwa mem -v 1 -t {threads} {input.ref[0]} {input.reads} 2> {log} | \
          samtools view --threads {threads} -bS - 2>> {log} | \
          samtools view --threads {threads} -uf 4 - 2>> {log} | \
          bamToFastq -i stdin -fq {output} >> {log} 2>&1
        """

localrules: preproc_link_filtered_reads
rule preproc_link_filtered_reads:
    input:
        lambda wildcards: expand(
            "Preprocessing/{{omic}}.{read}.{filtered}.fq",
            filtered=("trimmed." if wildcards.omic == "mg" else "trimmed.rna_filtered.") + ".".join([s + "_filtered" for s in get_filtering_genomes()]),
            read=["se"]
        )
    output:
        expand("Preprocessing/{{omic}}.{read}.preprocessed.fq", read=["r1", "r2", "se"])
    wildcard_constraints:
        omic="|".join(get_omics()),
    message:
        "Link filtered reads: {input}"
    run:
        for i,o in zip(input,output):
            shell("ln -fs $(echo {i} | cut -f 2 -d /) {o} && touch -h {o}")
