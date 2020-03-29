SAMTOOLS_MEM = str(round(float(BIGMEMCORE[:-1]) * 0.75 - 0.5)) + "G"

FILTERING_CMD = """
        TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "alignment_XXXXXX.bam")
        bwa mem -v 1 -t {threads} {input[1]} {input[0]} 2>> {log} | \
          samtools view --threads {threads} -bS - 2>> {log} | \
          samtools view --threads {threads} -uf 4 - 2>> {log} | \
          bamToFastq -i stdin -fq {output[0]} >> {log} 2>&1
        rm -rf $TMP_FILE
        """

SYMLINK_CMD = """
ln -fs $(echo {input[0]} | cut -f 2 -d /) {output[0]} && touch -h {output[0]}
ln -fs $(echo {input[1]} | cut -f 2 -d /) {output[1]} && touch -h {output[1]}
ln -fs $(echo {input[2]} | cut -f 2 -d /) {output[2]} && touch -h {output[2]}
"""


rule mg_filtering:
    input:
        'Preprocessing/mg.se.trimmed.fq',
        '{dir}/{ref}.fa'.format(dir=DBPATH + "/filtering", ref=FILTER),
        ancient(expand(
            "{p}/filtering/{filter}.{ext}", filter=FILTER,
            ext=['fa', 'fa.amb', 'fa.ann', 'fa.bwt', 'fa.pac', 'fa.sa'], p=DBPATH))
    output:
        'Preprocessing/mg.se.trimmed.%s_filtered.fq' % FILTER
    params:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: BIGCORENO
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/preprocessing_mg_filtering.log"
    message: "mg_filtering: Filtering mg reads against {FILTER}."
    shell:
        FILTERING_CMD


rule symlink_mg_preprocessed_files:
    input:
        'Preprocessing/mg.r1.fq',
        'Preprocessing/mg.r2.fq',
        'Preprocessing/mg.se.trimmed.%s_filtered.fq' % FILTER
    output:
        'Preprocessing/mg.r1.preprocessed.fq',
        'Preprocessing/mg.r2.preprocessed.fq',
        'Preprocessing/mg.se.preprocessed.fq'
    params:
        runtime = "1:00:00",
        mem = MEMCORE
    threads: 1
#    log: "logs/preprocessing_symlink_mg_preprocessed_files.log"
    message: "symlink_mg_preprocessed_files: Symlinking filtered mg reads."
    shell:
        SYMLINK_CMD
