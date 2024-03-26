import math

rule call_contig_depth:
    input:
        bam = "Assembly/{type}.reads.sorted.bam",
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        assembly_fai = "Assembly/%s.assembly.merged.fa.fai" % ASS,
        assembly_bed = "Assembly/%s.assembly.merged.fa.bed3" % ASS
    output:
        glob_dist = f"Stats/{{type}}/{ASS}.assembly.mosdepth.global.dist.txt",
        summary = f"Stats/{{type}}/{ASS}.assembly.mosdepth.summary.txt",
        region_dist = f"Stats/{{type}}/{ASS}.assembly.mosdepth.region.dist.txt",
        per_base = f"Stats/{{type}}/{ASS}.assembly.per-base.bed.gz",
        regions = f"Stats/{{type}}/{ASS}.assembly.regions.bed.gz",
        breadth = f"Stats/{{type}}/{ASS}.assembly.thresholds.bed.gz",
        read_count = f"Stats/{{type}}/{ASS}.assembly.contig_read_count.txt",
        contig_depth = f"Stats/{{type}}/{ASS}.assembly.contig_depth.txt",
        report = report(f"Stats/{{type}}/{ASS}.assembly.contig_flagstat.txt", category="Assembly")
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
        ass_dir = "Assembly"
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/analysis/call_contig_depth.{type}.log"
    message:
        "call_contig_depth: Getting data on assembly coverage with {wildcards.type} reads."
    shell:
        """
        ### Calculate coverage ###
        echo "Calculating contig base coverage ..." >> {log} 2>&1
        export MOSDEPTH_PRECISION=6
        mosdepth --fast-mode \
                 --by {input.assembly_bed} \
                 --threads 4 \
                 --thresholds 1,10,30 \
                 {params.stats_dir}/{params.ass_type}.assembly \
                 {input.bam} >> {log} 2>&1

        ### Calculate read counts ###
        echo "Calculating contig read counts ..." >> {log} 2>&1
        hts_nim_tools count-reads --threads 4 \
                                  {input.assembly_bed} \
                                  {input.bam} 2>> {log}\
          | sort -k1,1V > {output.read_count} 2>> {log} 

        ### Calculate average read coverage depth ###
        echo "Running BEDTools for average depth in each position ..." >> {log}
        export TMPDIR="{params.tmp_dir}"
        TMP_DEPTH=$(mktemp "depth_file_XXXXXX.txt")

        genomeCoverageBed -ibam {input.bam} 2>> {log} | grep -v "genome" > ${{TMP_DEPTH}} 2>> {log}

        echo "Calculating average contig coverage ..." >> {log}
        # This method of depth calculation was adapted and modified from the CONCOCT code
        perl {SRCDIR}/calcAvgCoverage.pl ${{TMP_DEPTH}} {input.assembly} > {output.contig_depth} 2>> {log}

        ### Clean up ###
        echo "Removing temporary files ..." >> {log}
        rm -f ${{TMP_DEPTH}} >> {log} 2>&1

        ### Calculate mappoing stats ###
        echo "Running samtools flagstat ..." >> {log}
        samtools flagstat {input.bam} 2>> {log} | cut -f1 -d ' ' > {output.report} 2>> {log}

        echo "Done." >> {log}
        """
