import math

rule call_gene_depth:
    input:
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        "Assembly/{type}.reads.sorted.bam",
        "Stats/%s/%s.assembly.length.txt" % (ASS,ASS)
    output:
        "Stats/{type}/annotation/{type}.gene_depth.hist",
        "Stats/{type}/annotation/{type}.gene_depth.avg",
        "Stats/{type}/annotation/{type}.gene.len"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        # Ensure at least 20gb of mem
        max(int(math.ceil(20 / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
    params:
        tmp_dir=TMPDIR
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/analysis_call_gene_depth.{type}.log"
    message:
        "call_gene_depth: Getting data on gene coverage with {wildcards.type} reads."
    shell:
        """
        export TMPDIR="{params.tmp_dir}"
        TMP_FILE=$(mktemp "annotation_cov_XXXXXX.bed")
        sortBed -g {input[2]} -i  {input[0]} | awk '$4 < $5' | coverageBed -hist -sorted -g {input[2]} -b {input[1]} -a "stdin" | grep -v "^all" > $TMP_FILE

        paste <(cat $TMP_FILE | cut -f9 | cut -f1 -d \";\" | sed -e \"s/ID=//g\") \
         <(cut -f10,11,12,13 $TMP_FILE) > {output[0]}
        rm $TMP_FILE

        # Record average gene coverage
        {SRCDIR}/calcAvgCoverage.pl {output[0]} > {output[1]}

        # Record gene length file
        cut -f 1,4 {output[0]} | uniq > {output[2]}
        """
