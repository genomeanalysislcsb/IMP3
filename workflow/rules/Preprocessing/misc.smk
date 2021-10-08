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

rule preproc_gzip_fastq_raw:
    input:
        "Preprocessing/{omic}.{reads}.fq"
    output:
        "Preprocessing/{omic}.{reads}.fq.gz"
    log:
        "logs/preproc_gzip_fastq_raw.{omic}.{reads}.log"
    benchmark:
        "benchmark/preproc_gzip_fastq_raw.{omic}.{reads}.tsv"
    wildcard_constraints:
        omic="|".join(get_omics()),
        reads="|".join(get_read_types_raw()),
    threads: 5
    resources:
        runtime="04:00:00",
        mem=get_mem_norm,
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "Archive: {input}"
    shell:
        "pigz -p {threads} -q {input}"

use rule preproc_gzip_fastq_raw as preproc_gzip_fastq_proc with:
    input:
        "Preprocessing/{omic}.{reads}.{stage}.fq"
    output:
        "Preprocessing/{omic}.{reads}.{stage}.fq.gz"
    log:
        "logs/preproc_gzip_fastq_proc.{omic}.{reads}.{stage}.log"
    benchmark:
        "benchmark/preproc_gzip_fastq_proc.{omic}.{reads}.{stage}.tsv"
    wildcard_constraints:
        omic="|".join(get_omics()),
        reads="|".join(set(get_read_types_trimmed() + get_read_types_trimmoed() + get_read_types_proc())),
