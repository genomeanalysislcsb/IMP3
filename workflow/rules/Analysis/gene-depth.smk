import math

rule call_gene_depth:
    input:
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        bam = "Assembly/{type}.reads.sorted.bam",
        assembly = f"Assembly/{ASS}.assembly.merged.fa",
        assembly_bed = f"Assembly/{ASS}.assembly.merged.fa.bed3"
    output:
        gene_bed = f"Analysis/annotation/{{type}}.annotation_CDS_RNA_hmms.bed",
        bed_hist = "Stats/{type}/annotation/{type}.gene_depth.hist",
        gene_depth = "Stats/{type}/annotation/{type}.gene_depth.avg",
        gene_len = "Stats/{type}/annotation/{type}.gene.len",
        glob_dist = f"Stats/{{type}}/annotation/{{type}}.mosdepth.global.dist.txt",
        summary = f"Stats/{{type}}/annotation/{{type}}.mosdepth.summary.txt",
        region_dist = f"Stats/{{type}}/annotation/{{type}}.mosdepth.region.dist.txt",
        per_base = f"Stats/{{type}}/annotation/{{type}}.per-base.bed.gz",
        regions = f"Stats/{{type}}/annotation/{{type}}.regions.bed.gz",
        breadth = f"Stats/{{type}}/annotation/{{type}}.thresholds.bed.gz",
        read_count = f"Stats/{{type}}/annotation/{{type}}.gene_read_count.txt",
        # gene_depth = f"Stats/{{type}}/annotation/{{type}}.gene_depth.txt"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        # Ensure at least 50gb of mem
        max(int(math.ceil(50 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 4)
    params:
        tmp_dir = TMPDIR,
        ass_type = ASS,
        ome = lambda wildcards: wildcards.type,
        stats_dir = lambda wildcards: f"Stats/{wildcards.type}/annotation/{wildcards.type}"
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/analysis/call_gene_depth.{type}.log"
    message:
        "call_gene_depth: Getting data on gene coverage with {wildcards.type} reads."
    shell:
        """
        ### Make feature bed ###
        echo "Creating gene bed file ..." >> {log} 2>&1
        # TODO: Change to working code!
        gxf2bed --input {input.gff} \
                --output {output.gene_bed}\
                --parent "ID" \
                --child "ID" \
                --feature "ID" \
                --threads {threads} >> {log} 2>&1
        
        ### Calculate coverage ###
        echo "Calculating gene base coverage ..." >> {log} 2>&1
        export MOSDEPTH_PRECISION=6
        mosdepth --fast-mode \
                 --by {output.gene_bed} \
                 --threads 4 \
                 --thresholds 1,10,30 \
                 {params.stats_dir}/{params.ome} \
                 {input.bam} >> {log} 2>&1
        
        ### Calculate read counts ###
        echo "Calculating gene read counts ..." >> {log} 2>&1
        hts_nim_tools count-reads --threads 4 \
                                  {output.gene_bed} \
                                  {input.bam} 2>> {log} \
          | sort -k1,1V > {output.read_count} 2>> {log} 

        ### Calculate average read coverage depth ###
        echo "Running BEDTools for average depth in each position ..." >> {log}
        # TODO: Make this work

        echo "Calculating average gene coverage ..." >> {log}
        # This method of depth calculation was adapted and modified from the CONCOCT code
        perl {SRCDIR}/calcAvgCoverage.pl XXXXXXXXXXX > {output.gene_depth} 2>> {log}

        echo "Done." >> {log}
        """


    # shell:
    #     """
    #     export TMPDIR="{params.tmp_dir}"
    #     TMP_FILE=$(mktemp "annotation_cov_XXXXXX.bed")

    #     echo "Calculating gene coverage ..." >> {log}
    #     bedtools sort -g {input.assembly_bed} \
    #                   -i {input.gff} 2>> {log} \
    #       | awk '$4 < $5' 2>> {log} \
    #       | bedtools coverage -a "stdin" \
    #                           -b {input.bam} \
    #                           -g <(cut -f 1,3 {input.assembly_bed}) \
    #                           -sorted \
    #                           -hist 2>> {log} \
    #       | grep -v "^all" > $TMP_FILE

    #     paste <(cut -f9 $TMP_FILE | cut -f1 -d \";\" | sed -e \"s/ID=//g\") \
    #           <(cut -f10,11,12,13 $TMP_FILE) > {output.bed_hist}
    #     rm $TMP_FILE

    #     # Record average gene coverage
    #     {SRCDIR}/calcAvgCoverage.pl {output.bed_hist} > {output.gene_depth} 2>> {log}

    #     # Record gene length file
    #     cut -f 1,4 {output.bed_hist} | uniq > {output.gene_len}
    #     """
