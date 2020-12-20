SAMTOOLS_MEM = str(round(float(BIGMEMCORE[:-1]) * 0.75 - 0.5)) + "G"

FILTERING_CMD = """
        TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "alignment_XXXXXX.bam")
        bwa mem -v 1 -t {threads} {input[1]} {input[0]} 2>> {log} | \
          samtools view --threads {threads} -bS - 2>> {log} | \
          samtools view --threads {threads} -uf 4 - 2>> {log} | \
          bamToFastq -i stdin -fq {output} >> {log} 2>&1
        rm -rf $TMP_FILE
        """

rule mg_filtering:
    input:
        filtering_mg_SE_input,
        filtering_filter
    output:
        'Preprocessing/mg.se.trimmed.{filterstep}.fq'
    wildcard_constraints:
        filterstep = "|".join([".".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: BIGCORENO
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/preprocessing_mg_filtering.{filterstep}.log"
    message: "mg_filtering: Filtering mg reads to get {wildcards.filterstep}."
    shell:
        FILTERING_CMD


rule symlink_mg_preprocessed_files:
    input:
        expand("Preprocessing/mg.{read}.fq",read=["r1","r2"]),
        expand("Preprocessing/mg.{read}.trimmed.{filtered}.fq",filtered=".".join([s + "_filtered" for s in FILTER]),read=["se"])
    output:
        expand('Preprocessing/mg.{read}.preprocessed.fq',read=["r1","r2","se"])
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    threads: 1
    message: "symlink_mg_preprocessed_files: Symlinking filtered mg reads."
    run:
        for i,o in zip(input,output):
             shell("ln -fs $(echo {i} | cut -f 2 -d /) {o} && touch -h {o}")


