rule reads_filtering:
    input:
        reads=filtering_reads_pe_input,
        ref=filtering_filter
    output:
        expand('Preprocessing/{{type}}.{read}.{{filterstep}}.fq', read=['r1', 'r2', 'se'])
    log:
        "logs/preprocessing_{type}_filtering.{filterstep}.log"
    wildcard_constraints:
        type = "|".join(TYPES),
        filterstep = "|".join(["trimmed(\.rna_filtered)?\." + "\.".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = "04:00:00",
        mem_mb = BIGMEMCORE_MB
    threads:
        getThreads(BIGCORENO, bigmem=True)
    params:
        mem = str(max(int(BIGMEMCORE_MB * 0.75), 1)) + 'M'
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "reads_filtering: Filtering {wildcards.type} reads to get {wildcards.filterstep}."
    shell:
        """
        TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "alignment_XXXXXX.bam")
        
        echo "Mapping reads." >> {log} 2>&1
        bwa mem -v 1 -t {threads} {input.ref[0]} {input.reads[0]} {input.reads[1]} 2>> {log} | \
         samtools view --threads {threads} -bS - > ${{TMP_FILE}} 2>> {log}
         
        echo "Removing mapped reads." >> {log} 2>&1
        samtools merge --threads {threads} -u - \
         <(samtools view --threads {threads} -u  -f 4 -F 264 ${{TMP_FILE}} 2>> {log}) \
         <(samtools view --threads {threads} -u -f 8 -F 260 ${{TMP_FILE}} 2>> {log}) \
         <(samtools view --threads {threads} -u -f 12 -F 256 ${{TMP_FILE}} 2>> {log}) 2>> {log} | \
         samtools view --threads {threads} -bF 0x800 - 2>> {log} | \
         samtools sort --threads {threads} -m {params.mem} -n - 2>> {log} | \
         bamToFastq -i stdin -fq {output[0]} -fq2 {output[1]} &>> {log}
        rm ${{TMP_FILE}}
        
        if [[ -s {input.reads[2]} ]]; then
            echo "Mapping reads." >> {log} 2>&1
            bwa mem -v 1 -t {threads} {input.ref[0]} {input.reads[2]} 2>> {log} | \
             samtools view --threads {threads} -bS - 2>> {log} | \
             samtools view --threads {threads} -uf 4 - 2>> {log} | \
             bamToFastq -i stdin -fq {output[2]} &>> {log}
        else
            echo "There are no singletons reads. {input.reads[2]} is empty. Generating empty {output[2]}" >> {log}
            touch {output[2]}
        fi
        """

localrules: symlink_preprocessed_reads_files
rule symlink_preprocessed_reads_files:
    input:
        lambda wildcards: expand(
            "Preprocessing/{{type}}.{read}.{filtered}.fq",
            filtered=('trimmed.' if wildcards.type == 'mg' else 'trimmed.rna_filtered.') + \
                ".".join([s + "_filtered" for s in FILTER]),
            read=['r1', 'r2', 'se']
        )
    output:
        expand('Preprocessing/{{type}}.{read}.preprocessed.fq', read=['r1', 'r2', 'se'])
    wildcard_constraints:
        type = "|".join(TYPES),
    threads: 1
    message:
        "symlink_preprocessed_reads_files: link to preprocessed {wildcards.type} reads"
    run:
        for i,o in zip(input,output):
            shell("ln -fs $(echo {i} | cut -f 2 -d /) {o} && touch -h {o}")

