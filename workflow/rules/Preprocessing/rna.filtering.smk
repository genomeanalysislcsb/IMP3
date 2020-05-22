FASTAFILES = expand("{path}/{files}.fasta", files=config["sortmerna"]["files"], path=DBPATH + "/sortmerna")

FASTAINDEXED = expand("{path}/idx/{files}", files=config["sortmerna"]["files"], path=DBPATH + "/sortmerna")

REF = ':'.join('%s,%s' % (a, b) for a, b in zip(FASTAFILES, FASTAINDEXED))

SORTME_MEM = int(MEMCORE[:-1]) * 1000

FILTER_RNA_SHELL = """
TMP_R12=$(mktemp --tmpdir={TMPDIR} -t "tmp_XXXXXX.R12")
merge-paired-reads.sh {input[0]} {input[1]} $TMP_R12.fq &> {log}
sortmerna --ref {params.ref} --reads $TMP_R12.fq --fastx --aligned $TMP_R12.rrna \
--other $TMP_R12.rrna_filt -a {threads} --paired_in -m {params.sortmernamem} --log &>> {log}
unmerge-paired-reads.sh $TMP_R12.rrna.fq {output[3]} {output[4]} &>> {log}
unmerge-paired-reads.sh $TMP_R12.rrna_filt.fq {output[0]} {output[1]} &>> {log}
sortmerna --ref {params.ref} --reads {input[2]} --fastx --aligned Preprocessing/mt.se.trimmed.rna \
--other Preprocessing/mt.se.trimmed.rna_filtered -a {threads} -m {params.sortmernamem} --log &>> {log}
cp $TMP_R12.rrna.log {output.stats_pe}
rm $TMP_R12*
"""

rule filter_rna:
    input:
        'Preprocessing/mt.r1.trimmed.fq',
        'Preprocessing/mt.r2.trimmed.fq',
        'Preprocessing/mt.se.trimmed.fq',
        expand(
            "{db}/sortmerna/idx/{files}.{ext}",
            files=config["sortmerna"]["files"],
            ext=['bursttrie_0.dat', 'kmer_0.dat', 'pos_0.dat', 'stats'],
            db=DBPATH
        )
    output:
        # not aligned
        norna_r1='Preprocessing/mt.r1.trimmed.rna_filtered.fq',
        norna_r2='Preprocessing/mt.r2.trimmed.rna_filtered.fq',
        norna_se='Preprocessing/mt.se.trimmed.rna_filtered.fq',
        # aligned
        rna_r1='Preprocessing/mt.r1.trimmed.rna.fq',
        rna_r2='Preprocessing/mt.r2.trimmed.rna.fq',
        rna_se='Preprocessing/mt.se.trimmed.rna.fq',
        # stats
        stats_pe='Preprocessing/mt.pe.trimmed.rna.log',
        stats_se='Preprocessing/mt.se.trimmed.rna.log',
    log:
        "logs/preprocessing_filter_rna.log"
    threads: 8
    params:
        ref = REF,
        sortmernamem = str(SORTME_MEM),
        runtime = "48:00:00",
        mem = MEMCORE
    conda:
        ENVDIR + "/IMP_preprocessing.yaml"
    message:
        "filter_rna: filtering rRNA from mt reads."
    shell:
        FILTER_RNA_SHELL
