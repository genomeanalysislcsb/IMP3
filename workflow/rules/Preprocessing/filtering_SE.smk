checkpoint reads_filtering_se:
    input:
        reads = filtering_reads_se_input,
        ref = filtering_filter
    output:
        temp('Preprocessing/{type}.se.{filterstep}.fq.gz')
    log:
        "logs/preprocessing/{type}_filtering.{filterstep}.log"
    wildcard_constraints:
        type = "|".join(TYPES),
        filterstep = "|".join(["trimmed(__rrna_proc)?__" + "__".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    params:
        mem = str(max(int(MAIN_BIGMEMCORE_MB * 0.75), 1)) + 'M',
        tmp_dir=TMPDIR
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "reads_filtering: Filtering {wildcards.type} reads to get {wildcards.filterstep}."
    shadow: "copy-minimal"
    shell:
        """
        bwa-mem2 mem -v 1 -K {BWA_BATCH} \
                     -t {threads} \
                     {input.ref[0]} \
                     {input.reads} 2> {log} | \
        samtools fastq -0 {output} \
                       -f 4 \
                       --threads 32 \
                       -n \
                       - 2>> {log}
        """

rule make_preprocessed_reads_files_se:
    input:
        # check = se_filter_checkpoint,
        fake_pe = expand("Preprocessing/{{type}}.{read}.fq", read=['r1', 'r2']),
        final_filter = lambda wildcards: expand("Preprocessing/{{type}}.{read}.{filtered}.fq.gz",
                                 filtered=('trimmed.' if wildcards.type == 'mg' \
                                           else 'trimmed__rrna_proc__') \
                                           + "__".join([s + "_filtered" for s in FILTER]),
                                 read=['se'])
    output:
        expand('Preprocessing/{{type}}.{read}.preprocessed.fq.gz', read=['se'])
    wildcard_constraints:
        type = "|".join(TYPES),
    threads:
        getThreads(4)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log:
        "logs/preprocessing/{type}_make_preprocessed_reads_files.log"
    message:
        "make_preprocessed_reads_files_se: Create final preprocessed {wildcards.type} reads file"
    shadow: "minimal"
    shell:
        """
        echo "Creating final preprocessed {wildcards.type} reads file ..." > {log}
        cp {input.final_filter} {output} >> {log} 2>&1
        echo "Done." >> {log}
        """
