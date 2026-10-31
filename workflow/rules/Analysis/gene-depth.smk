import math

ruleorder: index_merged_bam > call_gene_depth

rule call_gene_depth:
    input:
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        bam = "Assembly/{type}.reads.sorted.bam",
        bam_index = "Assembly/{type}.reads.sorted.bam.bai",
        assembly = f"Assembly/{ASS}.assembly.merged.fa",
        assembly_bed = f"Assembly/{ASS}.assembly.merged.fa.bed3"
    output:
        gene_bed = f"Analysis/annotation/{{type}}.annotation_CDS_RNA_hmms.bed.gz",
        gene_depth = "Stats/{type}/{type}.gene_depth.avg",
        gene_len = "Stats/{type}/{type}.gene.len",
        glob_dist = f"Stats/{{type}}/gene.mosdepth.global.dist.txt",
        summary = f"Stats/{{type}}/gene.mosdepth.summary.txt",
        region_dist = f"Stats/{{type}}/gene.mosdepth.region.dist.txt",
        per_base = f"Stats/{{type}}/gene.per-base.bed.gz",
        # regions = f"Stats/{{type}}/gene.regions.bed.gz",
        breadth = f"Stats/{{type}}/gene.thresholds.bed.gz",
        read_count = f"Stats/{{type}}/gene_read_count.txt"
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
        stats_dir = lambda wildcards: f"Stats/{wildcards.type}",
        bed_tmp = lambda wildcards: f"Analysis/annotation/{wildcards.type}.annotation_CDS_RNA_hmms.bed"
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/analysis/call_gene_depth.{type}.log"
    message:
        "call_gene_depth: Getting data on gene coverage with {wildcards.type} reads."
    shadow: "copy-minimal"
    shell:
        """
        # Ensure stats_dir is there
        mkdir -p {params.stats_dir} >> {log} 2>&1

        ### Make feature bed ###
        echo "Creating gene bed file ..." >> {log} 2>&1
        gxf2bed --input {input.gff} \
                --output {params.bed_tmp}\
                --parent "" \
                --child "" \
                --feature "ID" \
                --threads {threads} >> {log} 2>&1
        
        ### Calculate coverage ###
        echo "Calculating gene base coverage ..." >> {log} 2>&1
        export MOSDEPTH_PRECISION=6
        mosdepth --fast-mode \
                 --by {params.bed_tmp} \
                 --threads 4 \
                 --thresholds 1,10,30 \
                 {params.stats_dir}/gene \
                 {input.bam} >> {log} 2>&1
        
        ### Calculate read counts ###
        echo "Calculating gene read counts ..." >> {log} 2>&1
        hts_nim_tools count-reads --threads 4 \
                                  {params.bed_tmp} \
                                  {input.bam} 2>> {log} \
          | sort -k1,1V -k4,4V > {output.read_count} 2>> {log} 

        ### The average gene coverage is actually already in `gene.regions.bed.gz` so we just decompress it ###
        echo "Decompressing gene depth file ..." >> {log} 2>&1
        gunzip -c {params.stats_dir}/gene.regions.bed.gz > {output.gene_depth} 2>> {log}

        ### Writing gene length to file ###
        echo "Writing gene length to file ..." >> {log} 2>&1
        cut -f 11,4 {params.bed_tmp} | sed -e 's/,$//g' > {output.gene_len} 2>> {log}

        ### Compress bed file ###
        echo "Compressing bed file ..." >> {log} 2>&1
        gzip -f {params.bed_tmp} >> {log} 2>&1

        echo "Done." >> {log}
        """


    # shadow: "copy-minimal"
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
