rule reads_filtering:
    input:
        reads=filtering_reads_se_input,
        ref=filtering_filter
    output:
        temp('Preprocessing/{type}.se.{filterstep}.fq')
    log:
        "logs/preprocessing/{type}_filtering.{filterstep}.log"
    wildcard_constraints:
        type = "|".join(TYPES),
        filterstep = "|".join(["trimmed(\.rrna_filtered)?\." + "__".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "reads_filtering: Filtering {wildcards.type} reads to get {wildcards.filterstep}."
    shell:
        """
        bwa-mem2 mem -v 1 \
                     -t {threads} \
                     {input.ref[0]} \
                     {input.reads} 2> {log} | \
        samtools fastq -0 {output} \
                       -f 4 \
                       --threads 32 \
                       -n \
                       - 2>> {log}
        """

localrules: make_preprocessed_reads_files

rule make_preprocessed_reads_files:
    input:
        expand("Preprocessing/mg.{read}.fq", read=['r1', 'r2']),
        lambda wildcards: expand("Preprocessing/{{type}}.{read}.{filtered}.fq",
                                 filtered=('trimmed.' if wildcards.type == 'mg' \
                                           else 'trimmed__rrna_proc__') \
                                           + ".".join([s + "_filtered" for s in FILTER]),
                                 read=['se'])
    output:
        expand('Preprocessing/{{type}}.{read}.preprocessed.fq', read=['r1', 'r2', 'se'])
    wildcard_constraints:
        type = "|".join(TYPES),
    threads:
        1
    log:
        "logs/preprocessing/{type}_make_preprocessed_reads_files.log"
    message:
        "symlink_preprocessed_reads_files: link to preprocessed {wildcards.type} reads"
    run:
        for i, o in zip(input, output):
            shell("cp {i} {o} > {log} 2>&1")
