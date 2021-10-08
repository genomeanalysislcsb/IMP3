rule archive_annotation:
    input: # TODO change input
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        expand("Analysis/annotation/{omic}.{db}.counts.tsv", omic=get_omics(), db=get_hmms())
    output:
        "Analysis/annotation/intermediary.tar.gz"
    log:
        "logs/annotation_tar_annotation.log"
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_norm(1),
    params:
        folder="Analysis/annotation/intermediary/"
    message:
        "tar_annotation: Compressing intermediary steps of annotation."
    shell:
       "(tar cvzf {output} {params.folder} && rm -r {params.folder}) &> {log}"
