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
        tmp_dir=TMPDIR,
        method = filtering_method,
        deacon_abs = DEACON_ABS,
        deacon_rel = DEACON_REL,
        bwa_rule = BWA_RULE,
        script = os.path.join(SRCDIR, "host_filter_reads.py")
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "reads_filtering: Filtering {wildcards.type} reads to get {wildcards.filterstep}."
    shadow: "copy-minimal"
    shell:
        """
        if [[ "{params.method}" == "deacon_pass" ]]; then
            echo "Human reads were removed by Deacon in an earlier filter step: reads passed on." > {log}
            cp {input.reads} {output}
            exit 0
        fi
        if [[ "{params.method}" == "deacon" || "{params.method}" == "deacon_ref" ]]; then
            echo "Removing reads matching {input.ref[0]} with Deacon ($(deacon --version))." >> {log} 2>&1
            deacon filter --deplete -a {params.deacon_abs} -r {params.deacon_rel} -t {threads} \
                          -s {log}.single.json -o {output} {input.ref[0]} {input.reads} >> {log} 2>&1
            cat {log}.single.json >> {log}
            exit 0
        fi
        if [[ "{params.method}" == "bwa" ]]; then
            bwa-mem2 mem -v 1 -K {BWA_BATCH} -t {threads} {input.ref[0]} {input.reads} 2> {log} | \
                samtools view --threads {threads} -F 2304 - 2>> {log} | \
                python {params.script} {params.bwa_rule} --out {output} --stats {log}.single.json >> {log} 2>&1
            exit 0
        fi
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
        final_filter = lambda wildcards: expand("Preprocessing/{{type}}.{read}.{filtered}.fq.gz",
                                 filtered=('trimmed__' if wildcards.type == 'mg' \
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

rule deacon_index_build:
    '''Deacon index (default k 31, w 15) of a filter reference other than the host (PhiX, extra libraries), method deacon'''
    input:
        ancient(os.path.join(DBPATH, "filtering", "{ref}.fa"))
    output:
        os.path.join(DBPATH, "filtering", "deacon", "{ref}.k31w15.idx")
    wildcard_constraints:
        ref = "|".join(DEACON_REFS) if DEACON_REFS else "no_reference"
    log:
        "logs/preprocessing/deacon_index.{ref}.log"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        4
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "deacon_index_build: Deacon index of {wildcards.ref}."
    shell:
        """
        deacon index build -t {threads} {input} > {output}.$$.tmp 2> {log} && mv {output}.$$.tmp {output}
        """
