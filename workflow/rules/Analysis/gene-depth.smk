import math

def get_cores4min_mem(wildcards, attempt):
    min_mem = attempt * 50
    return max(int(math.ceil(min_mem / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)

rule call_gene_depth:
    input:
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        bam = "Assembly/{type}.reads.sorted.bam",
        assembly_bed = "Assembly/%s.assembly.merged.fa.bed3" % ASS
    output:
        bed_hist = "Stats/{type}/annotation/{type}.gene_depth.hist",
        gene_depth = "Stats/{type}/annotation/{type}.gene_depth.avg",
        gene_len = "Stats/{type}/annotation/{type}.gene.len"
    resources:
        runtime = SHORT_RUNTIME,
        # mem_mb = MAIN_BIGMEMCORE_MB
        mem_mb = HIGH_MEM_BIGMEMCORE_MB
    # threads:
    #     get_cores4min_mem
    # retries:
    #     6
    threads:
        getThreads(128, bigmem=True, high_mem_cluster=True)
    params:
        tmp_dir=TMPDIR
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/analysis/call_gene_depth.{type}.log"
    message:
        "call_gene_depth: Getting data on gene coverage with {wildcards.type} reads."
    shell:
        """
        export TMPDIR="{params.tmp_dir}"
        TMP_FILE=$(mktemp "annotation_cov_XXXXXX.bed")

        echo "Calculating gene coverage ..." >> {log}
        bedtools sort -g {input.assembly_bed} \
                      -i {input.gff} 2>> {log} \
          | awk '$4 < $5' 2>> {log} \
          | bedtools coverage -a "stdin" \
                              -b {input.bam} \
                              -g <(cut -f 1,3 {input.assembly_bed}) \
                              -sorted \
                              -hist 2>> {log} \
          | grep -v "^all" > $TMP_FILE

        paste <(cut -f9 $TMP_FILE | cut -f1 -d \";\" | sed -e \"s/ID=//g\") \
              <(cut -f10,11,12,13 $TMP_FILE) > {output.bed_hist}
        rm $TMP_FILE

        # Record average gene coverage
        {SRCDIR}/calcAvgCoverage.pl {output.bed_hist} > {output.gene_depth} 2>> {log}

        # Record gene length file
        cut -f 1,4 {output.bed_hist} | uniq > {output.gene_len}
        """
