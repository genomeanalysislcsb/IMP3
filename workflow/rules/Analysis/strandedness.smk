# featureCountsStranding: <type>: auto -> the strand of the library is detected from the reads mapped to the assembly:
# featureCounts on the CDS with -s 1 and -s 2; the larger count wins when it is >= 3 x the
# smaller (-s 1 or 2), otherwise the library is taken as unstranded (-s 0). Stats/<type>/<type>.strandedness.txt
# holds the chosen value for the featureCounts rules, .json the counts.
if "STRANDEDNESS_RULES" not in globals():
    STRANDEDNESS_RULES = True

    rule detect_strandedness:
        input:
            map = "Assembly/{type}.reads.sorted.bam",
            gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
        output:
            txt = "Stats/{type}/{type}.strandedness.txt",
            json = "Stats/{type}/{type}.strandedness.json"
        params:
            paired = lambda wildcards: "-p --countReadPairs" if ((wildcards.type == "mg" and len(MG) >= 2)
                                                                  or (wildcards.type == "mt" and len(MT) >= 2)) else "",
            min_ratio = 3
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            getThreads(8)
        conda:
            os.path.join(ENVDIR, "IMP_annotation.yaml")
        log:
            "logs/analysis/featurecount/strandedness.{type}.log"
        message:
            "detect_strandedness: strand of the {wildcards.type} library."
        shadow: "minimal"
        shell:
            """
            for s in 1 2; do
                featureCounts {params.paired} -t CDS -g ID -s $s -a {input.gff} -T {threads} -o fc.s$s.tsv {input.map} >> {log} 2>&1
            done
            s1=$(awk '$1 == "Assigned" {{print $2}}' fc.s1.tsv.summary)
            s2=$(awk '$1 == "Assigned" {{print $2}}' fc.s2.tsv.summary)
            choice=$(awk -v a="$s1" -v b="$s2" -v r={params.min_ratio} 'BEGIN {{ if (b >= r * a && b > 0) print 2; else if (a >= r * b && a > 0) print 1; else print 0 }}')
            echo $choice > {output.txt}
            printf '{{"type": "%s", "assigned_s1": %s, "assigned_s2": %s, "min_ratio": %s, "stranding": %s}}\\n' \\
                {wildcards.type} $s1 $s2 {params.min_ratio} $choice > {output.json}
            echo "strand: -s 1 $s1, -s 2 $s2 -> -s $choice" >> {log}
            """
