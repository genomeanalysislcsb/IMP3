rule reads_filtering:
    input:
        reads=filtering_reads_se_input,
        ref=filtering_filter
    output:
        'Preprocessing/{type}.se.{filterstep}.fq'
    log:
        "logs/preprocessing_{type}_filtering.{filterstep}.log"
    wildcard_constraints:
        type = "|".join(TYPES),
        filterstep = "|".join(["trimmed(\.rna_filtered)?\." + "\.".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = "08:00:00",
        mem_mb = BIGMEMCORE_MB
    threads: getThreads(BIGCORENO, bigmem=True)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "reads_filtering: Filtering {wildcards.type} reads to get {wildcards.filterstep}."
    shell:
        """
        bwa mem -v 1 -t {threads} {input.ref[0]} {input.reads} 2> {log} | \
          samtools view --threads {threads} -bS - 2>> {log} | \
          samtools view --threads {threads} -uf 4 - 2>> {log} | \
          bamToFastq -i stdin -fq {output} >> {log} 2>&1
        """

localrules: symlink_preprocessed_reads_files
rule symlink_preprocessed_reads_files:
    input:
        expand("Preprocessing/mg.{read}.fq", read=['r1', 'r2']),
        lambda wildcards: expand(
            "Preprocessing/{{type}}.{read}.{filtered}.fq",
            filtered=('trimmed.' if wildcards.type == 'mg' else 'trimmed.rna_filtered.') + \
                ".".join([s + "_filtered" for s in FILTER]),
            read=['se']
        )
    output:
        expand('Preprocessing/{{type}}.{read}.preprocessed.fq', read=['r1', 'r2', 'se'])
    wildcard_constraints:
        type = "|".join(TYPES),
    threads: 1
    message:
        "symlink_preprocessed_reads_files: link to preprocessed {wildcards.type} reads"
    run:
        for i,o in zip(input,output):
            shell("ln -fs $(echo {i} | cut -f 2 -d /) {o} && touch -h {o}")
