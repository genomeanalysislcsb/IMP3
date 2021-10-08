rule reads_count:
    input:
        lambda wildcards: get_read_files_raw(wildcards.omic) + get_read_files_trimmed(wildcards.omic, True) + get_read_files_proc(wildcards.omic)
    output:
        report("Stats/{omic}/{omic}.read_counts.txt", category="Preprocessing")
    log:
        "logs/preproc.reads_count.{omic}.log"
    benchmark:
        "benchmark/preproc.reads_count.{omic}.tsv"
    resources:
        runtime="00-02:00:00",
        mem=get_mem_norm
    message:
        "Count reads in {input}"
    shell:
        "for fname in {input}; do echo -e \"${{fname}}\\t$(($(wc -l < ${{fname}})/4))\" >> {output}; done > {output} 2> {log}"
