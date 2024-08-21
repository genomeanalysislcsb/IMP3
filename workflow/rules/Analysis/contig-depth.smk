import math

# ruleorder: index_bam > fasta2bed > index_assembly > call_contig_depth

# rule index_bam:
#     input:
#         bam = "Assembly/{type}.reads.sorted.bam"
#     output:
#         bam_index = "Assembly/{type}.reads.sorted.bam.bai"
#     resources:
#         runtime = SHORT_RUNTIME,
#         mem_mb = MAIN_BIGMEMCORE_MB
#     threads:
#         # Ensure at least 50gb of mem
#         max(int(math.ceil(50 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
#     conda:
#         ENVDIR + "/IMP_mapping.yaml"
#     log:
#         "logs/assembly/index_bam.{type}.log"
#     message:
#         "index_bam: Indexing {wildcards.type} bam."
#     shell:
#         """
#         echo "Indexing {input} ..." >> {log}
#         samtools index {input} >> {log} 2>&1

#         # Check output and print name
#         if [ -f {output} ]; then
#             echo "Indexing done. Index at: {output}." >> {log}
#         else
#             echo "Indexing failed." >> {log}
#             exit 1
#         fi

#         echo "Done." >> {log}
#         """

# rule index_assembly:
#     input:
#         f"Assembly/{ASS}.assembly.merged.fa"
#     output:
#         f"Assembly/{ASS}.assembly.merged.fa.fai"
#     resources:
#         runtime = SHORT_RUNTIME,
#         mem_mb = MAIN_MEMCORE_MB
#     threads:
#         4
#     conda:
#         os.path.join(ENVDIR, "IMP_mapping.yaml")
#     log:
#         "logs/assembly/index_assembly.log"
#     message:
#         "index_assembly: Indexing assembly."
#     shell:
#         """
#         echo "Indexing assembly ..." >> {log}
#         samtools faidx {input[0]} >> {log} 2>&1
#         echo "Done." >> {log}
#         """

# rule fasta2bed:
#     input:
#         assembly = f"Assembly/{ASS}.assembly.merged.fa",
#         assembly_index = f"Assembly/{ASS}.assembly.merged.fa.fai"
#     output:
#         f"Assembly/{ASS}.assembly.merged.fa.bed3"
#     resources:
#         runtime = VERY_SHORT_RUNTIME,
#         mem_mb = MAIN_MEMCORE_MB
#     threads:
#         4
#     log:
#         "logs/assembly/fasta2bed.log"
#     message:
#         "fasta2bed: Writing bed file for assembly."
#     shell:
#         """
#         echo "Writing bed3 file from assembly index ..." >> {log}
#         awk '{{print $1 \"\t0\t\" $2}}' {input.assembly_index} > {output} 2>> {log}
#         echo "Done." >> {log}
#         """

rule call_contig_depth:
    input:
        bam = "Assembly/{type}.reads.sorted.bam",
        bam_index = "Assembly/{type}.reads.sorted.bam.bai",
        assembly = f"Assembly/{ASS}.assembly.merged.fa",
        assembly_fai = f"Assembly/{ASS}.assembly.merged.fa.fai",
        assembly_bed = f"Assembly/{ASS}.assembly.merged.fa.bed3"
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
        runtime = MEDIUM_RUNTIME,
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
                 --threads 4 \
                 --by {input.assembly_bed} \
                 --thresholds 1,10,30 \
                 {params.stats_dir}/{params.ass_type}.assembly \
                 {input.bam} >> {log} 2>&1

        ### Calculate read counts ###
        echo "Calculating contig read counts ..." >> {log} 2>&1
        # TODO: The default SAM flag is 1796, which seems to occasionally lead to inconsistencies, 
        #       since a contig will have depth but 0 counts, maybe because its not a primary alignment
        hts_nim_tools count-reads --threads 4 \
                                  {input.assembly_bed} \
                                  {input.bam} 2>> {log} \
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
