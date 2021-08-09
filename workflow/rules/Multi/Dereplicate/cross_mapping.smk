rule cross_mapping:
    input:
        lambda wildcards: expand("{raw_directory}/Preprocessing/mg.{read}.preprocessed.fq",
                raw_directory=samples.path[wildcards.cross_sample],
                read=["r1","r2","se"]),
        lambda wildcards: ancient(expand("{raw_directory}/Assembly/{ass}.assembly.merged.{ext}",
                raw_directory=samples.path[wildcards.sample],
                ass=samples.assembly_type[wildcards.sample],
                ext=["fa.amb","fa.bwt","fa.pac","fa.sa","fa.ann","fa"]))
    output:
        "Dereplication/cross_mapping/{cross_sample}.on.{sample}.sorted.bam"
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: getThreads(BIGCORENO)
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/dereplicate_cross_mapping_{cross_sample}.on.{sample}.log"
    message: "cross_mapping: Mapping {wildcards.cross_sample} reads on {wildcards.sample} contigs."
    shell:
        """
        SAMHEADER="@RG\\tID:{wildcards.cross_sample}\\tSM:MG"
        PREFIX=Dereplication/cross_mapping/{wildcards.cross_sample}.on.{wildcards.sample}
        # merge paired and se
        samtools merge --threads {threads} -f $PREFIX.merged.bam \
         <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input[8]} {input[0]} {input[1]} 2>> {log}| \
         samtools view --threads {threads} -bS -) \
         <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input[8]} {input[2]} 2> {log}| \
         samtools view --threads {threads} -bS -) 2>> {log}
        # sort
        samtools sort --threads {threads} -m {SAMTOOLS_MEM} $PREFIX.merged.bam > $PREFIX.sorted.bam 2>> {log}
        rm $PREFIX.merged.bam
        """

rule call_contig_depth:
    input:
        "Dereplication/cross_mapping/{cross_sample}.on.{sample}.sorted.bam",
        lambda wildcards: expand("{raw_directory}/Assembly/{ass}.assembly.merged.{ext}",
                raw_directory=samples.path[wildcards.sample],
                ass=samples.assembly_type[wildcards.sample],
                ext=["fa","fa.fai","fa.bed3"])
    output:
        "Dereplication/cross_mapping/{sample}.contig_depth.{cross_sample}.txt",
        "Dereplication/cross_mapping/{sample}.contig_flagstat.{cross_sample}.txt"
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/analysis_call_contig_depth.{cross_sample}.on.{sample}.log"
    message: "call_contig_depth: Getting data on {wildcards.sample} assembly coverage with {wildcards.cross_sample} reads."
    shell:
        """
        TMP_DEPTH=$(mktemp --tmpdir={TMPDIR} -t "depth_file_XXXXXX.txt")
        genomeCoverageBed -ibam {input[0]} | grep -v "genome" > $TMP_DEPTH
        echo "Depth calculation done" >> {log}

        ## This method of depth calculation was adapted and modified from the CONCOCT code
        awk -v OFS='\t' 'BEGIN {{pc=""}}
        {{
        c=$1;
        if (c == pc) {{
        cov=cov+$2*$5;
        }} else {{
        print pc,cov;
        cov=$2*$5;
        pc=c}}
        }} END {{print pc,cov}}' $TMP_DEPTH | tail -n +2 > {output[0]}

        echo "Remove the temporary file" >> {log}
        rm $TMP_DEPTH
        echo "flagstat" >> {log}
        samtools flagstat {input[0]} 2>> {log} | cut -f1 -d ' ' > {output[1]}
        """


  
