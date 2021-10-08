# TODO: could create FastQC reports for all FASTQ files and get the number of reads from there
rule preproc_count_reads:
    input:
        lambda wildcards: get_read_files_raw(wildcards.omic) + get_read_files_trimmed(wildcards.omic, True) + get_read_files_proc(wildcards.omic)
    output:
        report("Stats/{omic}/{omic}.read_counts.txt", category="Preprocessing")
    log:
        "logs/preproc_count_reads.{omic}.log"
    benchmark:
        "benchmark/preproc_count_reads.{omic}.tsv"
    threads: 1
    resources:
        runtime="00-02:00:00",
        mem=get_mem_norm
    message:
        "Count reads: {input}"
    shell:
        "for fname in {input}; do echo -e \"${{fname}}\\t$(($(wc -l < ${{fname}})/4))\"; done > {output} 2> {log}"
