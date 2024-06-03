rule reads_filtering:
    input:
        reads = filtering_reads_pe_input,
        ref = filtering_filter
    output:
        temp(expand('Preprocessing/{{type}}.{read}.{{filterstep}}.fq', read=['r1', 'r2', 'se']))
    log:
        "logs/preprocessing/{type}_filtering.{filterstep}.log"
    wildcard_constraints:
        type = "|".join(TYPES),
        filterstep = "|".join(["trimmed(__rrna_proc)?__" + "__".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    params:
        mem = str(max(int(MAIN_BIGMEMCORE_MB * 0.75), 1)) + 'M',
        tmp_dir=TMPDIR
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "reads_filtering: Filtering {wildcards.type} reads to get {wildcards.filterstep}."
    shell:
        """
        export TMPDIR="{params.tmp_dir}"
        TMP_FILE=$(mktemp "alignment_XXXXXX.bam")
        
        echo "Mapping PE reads." >> {log} 2>&1

        bwa-mem2 mem -v 1 \
                     -t {threads} \
                     {input.ref[0]} \
                     {input.reads[0]} \
                     {input.reads[1]} 2>> {log} | \
         samtools view --threads {threads} \
                       --bam \
                       - > ${{TMP_FILE}} 2>> {log}
        
        echo "#####################################################################################" >> {log}
        
        echo "Removing mapped reads." >> {log} 2>&1
        samtools merge --threads {threads} \
                       -u \
                       -o - \
                      <(samtools view --threads {threads} \
                                      --uncompressed \
                                      --require-flags 4 \
                                      --exclude-flags 264 \
                                      ${{TMP_FILE}} 2>> {log}) \
                      <(samtools view --threads {threads} \
                                      --uncompressed \
                                      --require-flags 8 \
                                      --exclude-flags 260 \
                                      ${{TMP_FILE}} 2>> {log}) \
                      <(samtools view --threads {threads} \
                                      --uncompressed \
                                      --require-flags 12 \
                                      --exclude-flags 256 \
                                      ${{TMP_FILE}} 2>> {log}) 2>> {log} | \
        samtools view --threads {threads} \
                    --bam \
                    --exclude-flags 2048 \
                    - 2>> {log} | \
        samtools sort --threads {threads} \
                    -m {params.mem} \
                    -n \
                    - 2>> {log} | \
        samtools fastq -1 {output[0]} \
                       -2 {output[1]} \
                       -0 /dev/null \
                       -s /dev/null \
                       -f 4 \
                       --threads {threads} \
                       -n \
                       - >> {log} 2>&1
        rm ${{TMP_FILE}}

        echo "#####################################################################################" >> {log}

        if [[ -s {input.reads[2]} ]]; then
            echo "Mapping SE reads." >> {log} 2>&1

            bwa-mem2 mem -v 1 \
                        -t {threads} \
                        {input.ref[0]} \
                        {input.reads[2]} 2>> {log} | \
            samtools fastq -0 {output[2]} \
                           -f 4 \
                           --threads 32 \
                           -n \
                           - 2>> {log}
        else
            echo "There are no singletons reads. {input.reads[2]} is empty. Generating empty {output[2]}" >> {log}
            touch {output[2]}
        fi

        echo "#####################################################################################" >> {log}
        echo "Done." >> {log}
        """

localrules: make_preprocessed_reads_files

rule make_preprocessed_reads_files:
    input:
        lambda wildcards: expand("Preprocessing/{{type}}.{read}.{filtered}.fq",
                                 filtered=('trimmed__' if wildcards.type == 'mg' \
                                           else 'trimmed__rrna_proc__') \
                                           + "__".join([s + "_filtered" for s in FILTER]),
                                 read=['r1', 'r2', 'se'])
    output:
        expand('Preprocessing/{{type}}.{read}.preprocessed.fq', read=['r1', 'r2', 'se'])
    wildcard_constraints:
        type = "|".join(TYPES),
    threads:
        1
    log:
        "logs/preprocessing/{type}_make_preprocessed_reads_files.log"
    message:
        "make_preprocessed_reads_files: Create final preprocessed {wildcards.type} reads file"
    run:
        for i, o in zip(input,output):
            shell("cp {i} {o} > {log} 2>&1")

