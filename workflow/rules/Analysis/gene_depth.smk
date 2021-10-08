rule analysis_gene_depth_hist:
    input:
        gff="Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        bam="Assembly/{omic}.reads.sorted.bam",
        length="Stats/{asm}/{asm}.assembly.length.txt".format(asm=get_omics_assembly()),
    output:
        hist="Stats/{omic}/annotation/{omic}.gene_depth.hist",
        bed=temp("Stats/{omic}/annotation/{omic}.gene_depth.bed")
    log:
        "logs/analysis.gene_depth_hist.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics())
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_big(1)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "gene depth: hist: {input}"
    shell:
        "(sortBed -g {input.length} -i {input.gff} | awk '$4 < $5' | coverageBed -hist -sorted -g {input.length} -b {input.bam} -a stdin | grep -v '^all' > {output.bed} && "
        "paste <(cat {output.bed} | cut -f9 | cut -f1 -d ';' | sed -e 's/ID=//g') <(cut -f10,11,12,13 {output.bed}) > {output.hist}) 2> {log}"

rule analysis_gene_depth_avg:
    input:
        "Stats/{omic}/annotation/{omic}.gene_depth.hist"
    output:
        "Stats/{omic}/annotation/{omic}.gene_depth.avg",
    log:
        "logs/analysis.gene_depth_avg.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics())
    params:
        script=os.path.join(SRCDIR, "calcAvgCoverage.pl")
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_big(1)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "gene depth: avg: {input}"
    shell:
        "{params.script} {input} > {output} 2> {log}"

rule analysis_gene_depth_len:
    input:
        "Stats/{omic}/annotation/{omic}.gene_depth.hist"
    output:
        "Stats/{omic}/annotation/{omic}.gene.len"
    log:
        "logs/analysis.gene_depth_len.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics())
    threads: 1
    resources:
        runtime="02:00:00",
        mem=get_mem_norm(1)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "gene depth: len: {input}"
    shell:
        "(cut -f 1,4 {input} | uniq > {output}) 2> {log}"
