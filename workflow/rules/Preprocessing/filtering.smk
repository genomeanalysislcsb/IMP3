rule reads_filtering:
    input:
        reads=filtering_reads_pe_input,
        ref=filtering_filter
    output:
        expand("Preprocessing/{{omic}}.{read}.{{filterstep}}.fq", read=["r1", "r2", "se"])
    log:
        "logs/preprocessing_{omic}_filtering.{filterstep}.log"
    wildcard_constraints:
        omic="|".join(get_omics()),
        filterstep="|".join(["trimmed(\.rna_filtered)?\." + s for s in get_read_files_filtered_genomes_basenames(True)]),
    threads: get_cores_big()
    resources:
        runtime="08:00:00",
        mem=get_mem_big,
    params:
        mem=lambda wildcards, threads: str(round(get_mem_big(wildcards, None, threads, None, "Gb") * 0.75 - 0.5)) + "G" # TODO
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "reads_filtering: Filtering {wildcards.omic} reads to get {wildcards.filterstep}."
    shell:
        """
        TMP_FILE=$(mktemp --tmpdir=${{TMPDIR}} -t "alignment_XXXXXX.bam")
        
        bwa mem -v 1 -t {threads} {input.ref[0]} {input.reads[0]} {input.reads[1]} 2>> {log} | \
         samtools view --threads {threads} -bS - > ${{TMP_FILE}} 2>> {log}
        
        samtools merge --threads {threads} -u - \
         <(samtools view --threads {threads} -u  -f 4 -F 264 ${{TMP_FILE}} 2>> {log}) \
         <(samtools view --threads {threads} -u -f 8 -F 260 ${{TMP_FILE}} 2>> {log}) \
         <(samtools view --threads {threads} -u -f 12 -F 256 ${{TMP_FILE}} 2>> {log}) 2>> {log} | \
         samtools view --threads {threads} -bF 0x800 - 2>> {log} | \
         samtools sort --threads {threads} -m {params.mem} -n - 2>> {log} | \
         bamToFastq -i stdin -fq {output[0]} -fq2 {output[1]} &>> {log}
        rm ${{TMP_FILE}}
        
        if [[ -s {input.reads[2]} ]]; then
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
            "Preprocessing/{{omic}}.{read}.{filtered}.fq",
            filtered=("trimmed." if wildcards.omic == "mg" else "trimmed.rna_filtered.") + \
                ".".join([s + "_filtered" for s in get_filtering_genomes()]),
            read=["r1", "r2", "se"]
        )
    output:
        expand("Preprocessing/{{omic}}.{read}.preprocessed.fq", read=["r1", "r2", "se"])
    wildcard_constraints:
        omic="|".join(get_omics()),
    message:
        "symlink_preprocessed_reads_files: link to preprocessed {wildcards.omic} reads"
    run:
        for i,o in zip(input,output):
            shell("ln -fs $(echo {i} | cut -f 2 -d /) {o} && touch -h {o}")
