SAMTOOLS_MEM = str(round(float(BIGMEMCORE[:-1]) * 0.75 - 0.5)) + "G"

FILTERING_CMD = """
        TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "alignment_XXXXXX.bam")
        bwa mem -v 1 -t {threads} {input[3]} {input[0]} {input[1]} 2>> {log} | samtools view --threads {threads} -bS - > $TMP_FILE 2> {log}
        samtools merge --threads {threads} -u - \
         <(samtools view --threads {threads} -u  -f 4 -F 264 $TMP_FILE 2>> {log}) \
         <(samtools view --threads {threads} -u -f 8 -F 260 $TMP_FILE 2>> {log}) \
         <(samtools view --threads {threads} -u -f 12 -F 256 $TMP_FILE 2>> {log}) 2>> {log} | \
         samtools view --threads {threads} -bF 0x800 - 2> {log} | \
         samtools sort --threads {threads} -m {SAMTOOLS_MEM} -n - 2> {log} | \
         bamToFastq -i stdin -fq {output[0]} -fq2 {output[1]} >> {log} 2>&1
        if [[ -s {input[2]} ]]
         then
         bwa mem -v 1 -t {threads} {input[3]} {input[2]} 2>> {log} | \
          samtools view --threads {threads} -bS - 2>> {log} | \
          samtools view --threads {threads} -uf 4 - 2>> {log} | \
          bamToFastq -i stdin -fq {output[2]} >> {log} 2>&1
        else
         echo "{input[2]} is empty, skipping single end reference sequence filtering, but creating it anyway..." >> {log}
         touch {output[2]}
        fi
        rm -rf $TMP_FILE
        """

rule mg_filtering:
    input:
        filtering_mg_input,
        filtering_filter
    output:
        expand('Preprocessing/mg.{read}.trimmed.{{filterstep}}.fq',read=["r1","r2","se"])
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

rule mt_filtering:
    input:
        filtering_mt_input,
        filtering_filter
    output:
        expand('Preprocessing/mt.{read}.trimmed.rna_filtered.{{filterstep}}.fq',read=["r1","r2","se"])
    wildcard_constraints:
        filterstep = "|".join([".".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: BIGCORENO
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/preprocessing_mt_filtering.{filterstep}.log"
    message: "mt_filtering: Filtering mt reads to get {wildcards.filterstep}."
    shell:
        FILTERING_CMD

rule symlink_mt_preprocessed_files:
    input:
        expand("Preprocessing/mt.{read}.trimmed.rna_filtered.{filtered}.fq",filtered=".".join([s + "_filtered" for s in FILTER]),read=["r1","r2","se"])
    output:
        expand('Preprocessing/mt.{read}.preprocessed.fq',read=["r1","r2","se"])
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    threads: 1
    message: "symlink_mt_preprocessed_files: Symlinking filtered mt reads."
    run:
        for i,o in zip(input,output):
             shell("ln -fs $(echo {i} | cut -f 2 -d /) {o} && touch -h {o}")


rule symlink_mg_preprocessed_files:
    input:
        expand("Preprocessing/mg.{read}.trimmed.{filtered}.fq",filtered=".".join([s + "_filtered" for s in FILTER]),read=["r1","r2","se"])
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
