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
        runtime = "12:00:00",
        mem_mb = BIGMEMCORE_MB
    threads:
        # Ensure at least 5gb of mem
        max(int(math.ceil(5 / (int(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]]))),1)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/analysis_call_gene_depth.{type}.log"
    message:
        "call_gene_depth: Getting data on gene coverage with {wildcards.type} reads."
    shell:
        """
        TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "annotation_cov_XXXXXX.bed")
        sortBed -g {input[2]} -i  {input[0]} | awk '$4 < $5' | coverageBed -hist -sorted -g {input[2]} -b {input[1]} -a "stdin" | grep -v "^all" > $TMP_FILE

        paste <(cat $TMP_FILE | cut -f9 | cut -f1 -d \";\" | sed -e \"s/ID=//g\") \
         <(cut -f10,11,12,13 $TMP_FILE) > {output[0]}
        rm $TMP_FILE

        # Record average gene coverage
        {SRCDIR}/calcAvgCoverage.pl {output[0]} > {output[1]}

        # Record gene length file
        cut -f 1,4 {output[0]} | uniq > {output[2]}
        """
