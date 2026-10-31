rule filter_rrna_pe:
    input:
        r1 = "Preprocessing/mt.r1.trimmed.fq.gz",
        r2 = "Preprocessing/mt.r2.trimmed.fq.gz",
    output:
        filt_r1 = temp("Preprocessing/mt.r1.trimmed__rrna_filtered.fq.gz"),
        filt_r2 = temp("Preprocessing/mt.r2.trimmed__rrna_filtered.fq.gz"),
        rrna_r1 = "Preprocessing/mt.r1.trimmed.rrna.fq.gz",
        rrna_r2 = "Preprocessing/mt.r2.trimmed.rrna.fq.gz",
        # otu = "Preprocessing/otu_map.txt", # TODO: Check if this is needed
    log:
        "logs/preprocessing/filter_rrna_pe.log"
    threads:
        getThreads(128)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "filter_rna: filtering rRNA from mt PE reads."
    # shadow: "minimal"
    shell:
        """
        ### Get average read length ###
        echo "Calculating average read length ..." > {log}
        r_mean_l=$( seqkit stats -T {input.r1} 2>> {log} | cut -f 7 2>> {log} | tail -n +2 2>> {log} | awk '{{printf("%d\\n", $1 + 0.5)}}' 2>> {log} )

        ### Run RiboDetector ###
        echo "Running RiboDetector ..." >> {log}

        ribodetector_cpu --threads {threads} \
                         --len ${{r_mean_l}} \
                         --input {input.r1} {input.r2} \
                         --ensure none \
                         --chunk_size 256 \
                         --output {output.filt_r1} {output.filt_r2} \
                         --rrna {output.rrna_r1} {output.rrna_r2} >> {log} 2>&1

        echo "Done." >> {log}
        """

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
        r_mean_l=$( seqkit stats -T {input.se} 2>> {log} | cut -f 7 2>> {log} | tail -n +2 2>> {log} | awk '{{printf("%d\\n", $1 + 0.5)}}' 2>> {log} )

        ### Run RiboDetector ###
        echo "Running RiboDetector ..." >> {log}

        ribodetector_cpu --threads {threads} \
                         --len ${{r_mean_l}} \
                         --input {input.se} \
                         --ensure none \
                         --chunk_size 256 \
                         --output {output.filt_se} \
                         --rrna {output.rrna_se} >> {log} 2>&1
        
        echo "Done." >> {log}
        """

# If KEEP_RRNAs is set to True, then the rRNA-filtered reads are concatenated with rRNA reads
if KEEP_RRNAs:
    rule copy_rrna_proc_out:
        input:
            r1_filt = rules.filter_rrna_pe.output.filt_r1,
            r2_filt = rules.filter_rrna_pe.output.filt_r2,
            se_filt = rules.filter_rrna_se.output.filt_se,
            r1_rrna = rules.filter_rrna_pe.output.rrna_r1,
            r2_rrna = rules.filter_rrna_pe.output.rrna_r2,
            se_rrna = rules.filter_rrna_se.output.rrna_se
        output:
            r1 = temp("Preprocessing/mt.r1.trimmed__rrna_proc.fq.gz"),
            r2 = temp("Preprocessing/mt.r2.trimmed__rrna_proc.fq.gz"),
            se = temp("Preprocessing/mt.se.trimmed__rrna_proc.fq.gz")
        log:
            "logs/preprocessing/copy_rrna_proc_out.log"
        threads:
            getThreads(4)
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        message:
            "filter_rna: concatenating rRNA reads with other reads."
        shell:
            """
            echo "Concatenating rRNA-filtered reads with rRNA reads ..." > {log}
            cat {input.r1_filt} {input.r1_rrna} > {output.r1} 2>> {log}
            cat {input.r2_filt} {input.r2_rrna} > {output.r2} 2>> {log}
            cat {input.se_filt} {input.se_rrna} > {output.se} 2>> {log}

            echo "Done." >> {log}
            """
else:
    rule copy_rrna_proc_out:
        input:
            r1 = rules.filter_rrna_pe.output.filt_r1,
            r2 = rules.filter_rrna_pe.output.filt_r2,
            se = rules.filter_rrna_se.output.filt_se,
            r1_rrna = rules.filter_rrna_pe.output.rrna_r1,
            r2_rrna = rules.filter_rrna_pe.output.rrna_r2,
            se_rrna = rules.filter_rrna_se.output.rrna_se
        output:
            r1 = temp("Preprocessing/mt.r1.trimmed__rrna_proc.fq.gz"),
            r2 = temp("Preprocessing/mt.r2.trimmed__rrna_proc.fq.gz"),
            se = temp("Preprocessing/mt.se.trimmed__rrna_proc.fq.gz")
        log:
            "logs/preprocessing/copy_rrna_proc_out.log"
        threads:
            getThreads(4)
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        message:
            "filter_rna: copying rRNA-filtered reads."
        shell:
            """
            echo "Only keeping rRNA-filtered reads ..." > {log}
            cp {input.r1} {output.r1} 2>> {log}
            cp {input.r2} {output.r2} 2>> {log}
            cp {input.se} {output.se} 2>> {log}

            echo "Done." >> {log}
            """
