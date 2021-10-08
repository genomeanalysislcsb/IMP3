# Archive data from the "Preprocessing" step

rule archive_reads_raw:
    input:
        "Preprocessing/{omic}.{reads}.fq"
    output:
        "Preprocessing/{omic}.{reads}.fq.gz"
    log:
        "logs/archive_reads_raw.{omic}.{reads}.log"
    benchmark:
        "benchmark/archive_reads_raw.{omic}.{reads}.tsv"
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

use rule archive_reads_raw as archive_reads_proc with:
    input:
        "Preprocessing/{omic}.{reads}.{stage}.fq"
    output:
        "Preprocessing/{omic}.{reads}.{stage}.fq.gz"
    log:
        "logs/archive_reads_proc.{omic}.{reads}.{stage}.log"
    benchmark:
        "benchmark/archive_reads_proc.{omic}.{reads}.{stage}.tsv"
    wildcard_constraints:
        omic="|".join(get_omics()),
        reads="|".join(set(get_read_types_trimmed() + get_read_types_trimmoed() + get_read_types_proc())),
