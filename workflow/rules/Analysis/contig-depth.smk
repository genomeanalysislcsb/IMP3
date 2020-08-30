rule call_contig_depth:
    input:
        "Assembly/{type}.reads.sorted.bam",
        "Assembly/%s.assembly.merged.fa" % ASS,
        "Assembly/%s.assembly.merged.fa.fai" % ASS,
        "Assembly/%s.assembly.merged.fa.bed3" % ASS
    output:
        "Stats/{type}/%s.assembly.contig_coverage.txt" % ASS,
        "Stats/{type}/%s.assembly.contig_depth.txt" % ASS,
        report("Stats/{type}/%s.assembly.contig_flagstat.txt" % ASS,category="Assembly")
    params:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/analysis_call_contig_depth.{type}.log"
    message: "call_contig_depth: Getting data on assembly coverage with {wildcards.type} reads."
    shell:
        """
        coverageBed -b {input[0]} -a {input[3]} -sorted > {output[0]} 2>> {log}
        echo "Coverage calculation done" >> {log}
        echo "Running BEDTools for average depth in each position" >> {log}
        TMP_DEPTH=$(mktemp --tmpdir={TMPDIR} -t "depth_file_XXXXXX.txt")
        genomeCoverageBed -ibam {input[0]} | grep -v "genome" > $TMP_DEPTH
        echo "Depth calculation done" >> {log}

        ## This method of depth calculation was adapted and modified from the CONCOCT code
	{SRCDIR}/calcAvgCoverage.pl $TMP_DEPTH {input[1]} > {output[1]}	

        echo "Remove the temporary file" >> {log}
        rm $TMP_DEPTH
        echo "flagstat" >> {log}
        samtools flagstat {input[0]} 2>> {log} | cut -f1 -d ' ' > {output[2]}
        """
