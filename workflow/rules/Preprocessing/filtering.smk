rule reads_filtering:
    input:
        reads = filtering_reads_pe_input,
        ref = filtering_filter
    output:
        temp(expand('Preprocessing/{{type}}.{read}.{{filterstep}}.fq.gz', read=['r1', 'r2', 'se']))
    log:
        "logs/preprocessing/{type}_filtering.{filterstep}.log"
    wildcard_constraints:
        type = "|".join(TYPES),
        filterstep = "|".join(["trimmed(__rrna_proc)?__" + "__".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    params:
        mem = str(max(int(MAIN_BIGMEMCORE_MB * 0.75), 1)) + 'M',
        tmp_dir=TMPDIR,
        method = filtering_method,
        deacon_abs = DEACON_ABS,
        deacon_rel = DEACON_REL,
        bwa_rule = BWA_RULE,
        script = os.path.join(SRCDIR, "host_filter_reads.py")
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "reads_filtering: Filtering {wildcards.type} reads to get {wildcards.filterstep}."
    shadow: "copy-minimal"
    shell:
        """
        export TMPDIR="{params.tmp_dir}"
        if [[ "{params.method}" == "deacon_pass" ]]; then
            echo "Human reads were removed by Deacon in an earlier filter step: reads passed on." >> {log} 2>&1
            cp {input.reads[0]} {output[0]} && cp {input.reads[1]} {output[1]} && cp {input.reads[2]} {output[2]}
            exit 0
        fi
        if [[ "{params.method}" == "deacon" || "{params.method}" == "deacon_ref" ]]; then
            echo "Removing reads matching {input.ref[0]} with Deacon ($(deacon --version))." >> {log} 2>&1
            deacon filter --deplete -a {params.deacon_abs} -r {params.deacon_rel} -t {threads} \
                          -s {log}.paired.json -o {output[0]} -O {output[1]} \
                          {input.ref[0]} {input.reads[0]} {input.reads[1]} >> {log} 2>&1
            cat {log}.paired.json >> {log}
            if [[ -s {input.reads[2]} ]]; then
                deacon filter --deplete -a {params.deacon_abs} -r {params.deacon_rel} -t {threads} \
                              -s {log}.single.json -o {output[2]} {input.ref[0]} {input.reads[2]} >> {log} 2>&1
                cat {log}.single.json >> {log}
            else
                echo "There are no singletons reads. {input.reads[2]} is empty. Generating empty {output[2]}" >> {log}
                touch {output[2]}
            fi
            echo "Done." >> {log}
            exit 0
        fi
        if [[ "{params.method}" == "bwa" ]]; then
            echo "Mapping PE reads; a pair is removed when a mate aligns over most of its length ({params.bwa_rule})." >> {log}
            bwa-mem2 mem -v 1 -K {BWA_BATCH} -t {threads} {input.ref[0]} {input.reads[0]} {input.reads[1]} 2>> {log} | \
                samtools view --threads {threads} -F 2304 - 2>> {log} | \
                python {params.script} --paired {params.bwa_rule} --out1 {output[0]} --out2 {output[1]} \
                    --stats {log}.paired.json >> {log} 2>&1
            if [[ -s {input.reads[2]} ]]; then
                bwa-mem2 mem -v 1 -K {BWA_BATCH} -t {threads} {input.ref[0]} {input.reads[2]} 2>> {log} | \
                    samtools view --threads {threads} -F 2304 - 2>> {log} | \
                    python {params.script} {params.bwa_rule} --out {output[2]} --stats {log}.single.json >> {log} 2>&1
            else
                echo "There are no singletons reads. {input.reads[2]} is empty. Generating empty {output[2]}" >> {log}
                touch {output[2]}
            fi
            echo "Done." >> {log}
            exit 0
        fi
        # bwa_any (IMP3 <= 3.0.3): a pair is removed when bwa-mem2 reports either mate as mapped
        TMP_FILE=$(mktemp "alignment_XXXXXX.bam")
        
        echo "Mapping PE reads." >> {log} 2>&1

        bwa-mem2 mem -v 1 -K {BWA_BATCH} \
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

            bwa-mem2 mem -v 1 -K {BWA_BATCH} \
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

rule make_preprocessed_reads_files:
    input:
        lambda wildcards: expand("Preprocessing/{{type}}.{read}.{filtered}.fq.gz",
                                 filtered=('trimmed__' if wildcards.type == 'mg' \
                                           else 'trimmed__rrna_proc__') \
                                           + "__".join([s + "_filtered" for s in FILTER]),
                                 read=['r1', 'r2', 'se'])
    output:
        expand('Preprocessing/{{type}}.{read}.preprocessed.fq.gz', read=['r1', 'r2', 'se'])
    wildcard_constraints:
        type = "|".join(TYPES),
    threads:
        getThreads(8)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log:
        "logs/preprocessing/{type}_make_preprocessed_reads_files.log"
    message:
        "make_preprocessed_reads_files: Create final preprocessed {wildcards.type} reads file"
    run:
        for i, o in zip(input,output):
            shell("cp {i} {o} > {log} 2>&1")

rule deacon_index_build:
    '''Deacon index (default k 31, w 15) of a filter reference other than the host (PhiX, extra libraries), method deacon'''
    input:
        ancient(os.path.join(DBPATH, "filtering", "{ref}.fa"))
    output:
        os.path.join(DBPATH, "filtering", "deacon", "{ref}.k31w15.idx")
    wildcard_constraints:
        ref = "|".join(DEACON_REFS) if DEACON_REFS else "no_reference"
    log:
        "logs/preprocessing/deacon_index.{ref}.log"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        4
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "deacon_index_build: Deacon index of {wildcards.ref}."
    shell:
        """
        deacon index build -t {threads} {input} > {output}.$$.tmp 2> {log} && mv {output}.$$.tmp {output}
        """
