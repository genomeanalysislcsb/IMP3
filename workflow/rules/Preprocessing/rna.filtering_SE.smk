rule filter_rrna_se:
    input:
        se = "Preprocessing/mt.se.trimmed.fq.gz",
    output:
        filt_se = temp("Preprocessing/mt.se.trimmed__rrna_filtered.fq.gz"),
        rrna_se = "Preprocessing/mt.se.trimmed.rrna.fq.gz",
    log:
        "logs/preprocessing/filter_rrna_se.log"
    threads:
        getThreads(128)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "filter_rna: filtering rRNA from mt SE reads."
    shadow: "copy-minimal"
    shell:
        """
        ### Get average read length ###
        echo "Calculating average read length ..." > {log}
        r_mean_l=$( seqkit stats -T {input.se} 2>> {log} | cut -f 7 | tail -n +2 | awk '{{printf("%d\\n", $1 + 0.5)}}' )

        ### Run RiboDetector ###
        echo "Running RiboDetector ..." >> {log}
        ribodetector_cpu --threads {threads} \\
                         --len ${{r_mean_l}} \\
                         --input {input.se} \\
                         --ensure none \\
                         --chunk_size 256 \\
                         --output {output.filt_se} \\
                         --rrna {output.rrna_se} >> {log} 2>&1

        echo "Done." >> {log}
        """

# keep_rrnas: true -> the rRNA reads are added back to the rRNA-filtered reads
rule copy_rrna_proc_out:
    input:
        se = rules.filter_rrna_se.output.filt_se,
        se_rrna = rules.filter_rrna_se.output.rrna_se
    output:
        se = temp("Preprocessing/mt.se.trimmed__rrna_proc.fq.gz")
    params:
        keep = "yes" if KEEP_RRNAs else "no"
    log:
        "logs/preprocessing/copy_rrna_proc_out.log"
    threads:
        getThreads(4)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    message:
        "filter_rna: rRNA-processed mt SE reads."
    shell:
        """
        if [[ "{params.keep}" == "yes" ]]; then
            echo "Concatenating rRNA-filtered reads with rRNA reads ..." > {log}
            cat {input.se} {input.se_rrna} > {output.se} 2>> {log}
        else
            echo "Only copying rRNA-filtered reads ..." > {log}
            cp {input.se} {output.se} 2>> {log}
        fi
        echo "Done." >> {log}
        """
