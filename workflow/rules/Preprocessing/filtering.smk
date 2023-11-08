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
        filterstep = "|".join(["trimmed(__rna_filtered)?__" + "__".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = "04:00:00",
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
        
        echo "Mapping reads." >> {log} 2>&1

        tree {params.tmp_dir} >> {log}
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
         bamToFastq -i stdin \
                    -fq {output[0]} \
                    -fq2 {output[1]} &>> {log}
        rm ${{TMP_FILE}}

        echo "#####################################################################################" >> {log}

        if [[ -s {input.reads[2]} ]]; then
            echo "Mapping SE reads." >> {log} 2>&1

            # atm piping directly from bwa --> samtools view --> bamToFastq fails for this part specifically??
            TMP_FILE_SE=$(mktemp "al_se_XXXXX.sam")
            TMP_FILE_SE2=$(mktemp "al_se_XXXXX.bam")
            
            bwa-mem2 mem -v 1 \
                         -t {threads} \
                         {input.ref[0]} \
                         {input.reads[2]} > ${{TMP_FILE_SE}} 2>> {log}
            
            samtools view --threads {threads} \
                          --bam \
                          --uncompressed \
                          --require-flags 4 \
                          ${{TMP_FILE_SE}} > ${{TMP_FILE_SE2}} 2>> {log}
            
            bamToFastq -i ${{TMP_FILE_SE2}} \
                       -fq {output[2]} &>> {log}

            rm ${{TMP_FILE_SE}} ${{TMP_FILE_SE2}}
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
            filtered=('trimmed__' if wildcards.type == 'mg' else 'trimmed__rna_filtered__') + \
                "__".join([s + "_filtered" for s in FILTER]),
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

