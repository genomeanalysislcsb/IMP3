rule analysis_archive:
    input:
        "status/analysis.ready"
    output:
        "Analysis/annotation/intermediary.tar.gz"
    log:
        "logs/analysis_archive.log"
    benchmark:
        "benchmark/analysis_archive.tsv"
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_norm,
    params:
        folder="Analysis/annotation/intermediary/"
    message:
        "Archive intermediate analysis output files"
    shell:
       "(tar cvzf {output} {params.folder} && rm -r {params.folder}) &> {log}"
