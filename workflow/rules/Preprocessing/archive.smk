# Archive data from the "Preprocessing" step

rule archive_reads_raw:
    input:
        "Preprocessing/{omic}.{reads}.fq"
    output:
        "Preprocessing/{omic}.{reads}.fq.gz"
    log:
        "logs/preproc.archive_reads_raw.{omic}.{reads}.log"
    benchmark:
        "benchmark/preproc.archive_reads_raw.{omic}.{reads}.tsv"
    wildcard_constraints:
        omic = "|".join(get_omics()),
        reads = "|".join(get_read_types_raw("mg")),
    threads: 5
    resources:
        runtime = "00-04:00:00",
        mem = get_mem_norm(5),
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "Gzip {input}"
    shell:
        "pigz -p {threads} -q {input}"

use rule archive_reads_raw as archive_reads_proc with:
    input:
        "Preprocessing/{omic}.{reads}.{stage}.fq"
    output:
        "Preprocessing/{omic}.{reads}.{stage}.fq.gz"
    log:
        "logs/preproc.archive_reads_proc.{omic}.{reads}.{stage}.log"
    benchmark:
        "benchmark/preproc.archive_reads_proc.{omic}.{reads}.{stage}.tsv"
    wildcard_constraints:
        omic = "|".join(get_omics()),
        reads = "|".join(set(get_read_types_trimmed("mg") + get_read_types_trimmoed("mg") + get_read_types_proc("mg"))),
