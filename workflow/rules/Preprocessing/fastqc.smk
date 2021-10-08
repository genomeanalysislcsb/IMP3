rule preproc_fastqc_raw:
    input:
        "Preprocessing/{omic}.{reads}.fq"
    output:
        expand("Stats/{{omic}}/{{omic}}.{{reads}}_fastqc.{ext}", ext=["zip", "html"])
    log:
        "logs/preproc_fastqc_raw.{omic}.{reads}.log"
    benchmark:
        "benchmark/preproc_fastqc_raw.{omic}.{reads}.tsv"
    wildcard_constraints:
        omic="|".join(get_omics()),
        reads="|".join(["r1", "r2", "se"]),
    threads: 5
    resources:
        runtime="01:00:00",
        mem=get_mem_norm,
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "FastQC: {input}"
    shell:
        "(if [[ -s {input} ]]; then fastqc -o $(dirname {output[0]}) -f fastq {input} -t {threads}; else echo \"Empty input: {input}\"; touch {output}; fi) &> {log}"

use rule preproc_fastqc_raw as preproc_fastqc_proc with:
    input:
        "Preprocessing/{omic}.{reads}.preprocessed.fq"
    output:
        expand("Stats/{{omic}}/{{omic}}.{{reads}}.preprocessed_fastqc.{ext}", ext=["zip", "html"])
    log:
        "logs/preproc_fastqc_proc.{omic}.{reads}.log"
    benchmark:
        "benchmark/preproc_fastqc_proc.{omic}.{reads}.tsv"
