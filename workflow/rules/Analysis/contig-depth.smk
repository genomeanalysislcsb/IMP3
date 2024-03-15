import math

def get_cores4min_mem(wildcards, attempt):
    min_mem = attempt * 50
    return max(int(math.ceil(min_mem / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)

rule call_contig_depth:
    input:
        bam = "Assembly/{type}.reads.sorted.bam",
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        assembly_fai = "Assembly/%s.assembly.merged.fa.fai" % ASS,
        assembly_bed = "Assembly/%s.assembly.merged.fa.bed3" % ASS
    output:
        contig_cov = "Stats/{type}/%s.assembly.contig_coverage.txt" % ASS,
        contig_depth = "Stats/{type}/%s.assembly.contig_depth.txt" % ASS,
        report = report("Stats/{type}/%s.assembly.contig_flagstat.txt" % ASS, category="Assembly")
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
        "logs/analysis/call_contig_depth.{type}.log"
    message:veggy
        "call_contig_depth: Getting data on assembly coverage with {wildcards.type} reads."
    shell:
        """
        samtools idxstats {input.bam} 2>> {log} | cut -f 1-2 > Assembly/bedtools_contig_cov_contigs.txt
        # Then sort the BED file according to this genome file
        bedtools sort -faidx Assembly/bedtools_contig_cov_contigs.txt -i {input.assembly_bed} > Assembly/bedtools_contig_cov_features-sorted.bed 2>> {log}
        # Finally use coverageBed with "-sorted", sorted BED file, and genome file
        echo "Calculating contig coverage ..." >> {log}
        bedtools coverage -a Assembly/bedtools_contig_cov_features-sorted.bed \
                          -b {input.bam} \
                          -g Assembly/bedtools_contig_cov_contigs.txt \
                          -sorted > {output.contig_cov} 2>> {log}
        echo "Coverage calculation done." >> {log}

        echo "Running BEDTools for average depth in each position" >> {log}
        export TMPDIR="{params.tmp_dir}"
        TMP_DEPTH=$(mktemp "depth_file_XXXXXX.txt")

        genomeCoverageBed -ibam {input.bam} 2>> {log} | grep -v "genome" > $TMP_DEPTH
        echo "Depth calculation done" >> {log}

        ## This method of depth calculation was adapted and modified from the CONCOCT code
        perl {SRCDIR}/calcAvgCoverage.pl $TMP_DEPTH {input.assembly} > {output.contig_depth} 2>> {log}

        echo "Remove the temporary files" >> {log}
        rm -f $TMP_DEPTH Assembly/bedtools_contig_cov_contigs.txt Assembly/bedtools_contig_cov_features-sorted.bed

        echo "flagstat" >> {log}
        samtools flagstat {input.bam} 2>> {log} | cut -f1 -d ' ' > {output.report}
        """
