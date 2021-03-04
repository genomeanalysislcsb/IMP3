FASTAFILES = expand("{path}/{files}.fasta", files=config["rnafiltering"]["files"], path=DBPATH + "/sortmerna")

FASTAINDEXED = expand("{path}/idx/{files}", files=config["rnafiltering"]["files"], path=DBPATH + "/sortmerna")

REF = ':'.join('%s,%s' % (a, b) for a, b in zip(FASTAFILES, FASTAINDEXED))

SORTME_MEM = int(MEMCORE[:-1]) * 1000

FILTER_RNA_SHELL = """
  TMP_R12=$(mktemp --tmpdir={TMPDIR} -t "tmp_XXXXXX.R12")
  merge-paired-reads.sh {input[0]} {input[1]} $TMP_R12.fq > {log} 2>&1
  sortmerna --ref {params.ref} --reads $TMP_R12.fq --fastx --aligned $TMP_R12.rrna \
   --other $TMP_R12.rrna_filt -a {threads} --paired_in -m {params.sortmernamem} >> {log} 2>&1
  unmerge-paired-reads.sh $TMP_R12.rrna.fq {output[3]} {output[4]} >> {log} 2>&1
  unmerge-paired-reads.sh $TMP_R12.rrna_filt.fq {output[0]} {output[1]} >> {log} 2>&1
  sortmerna --ref {params.ref} --reads {input[2]} --fastx --aligned Preprocessing/mt.se.trimmed.rna \
   --other Preprocessing/mt.se.trimmed.rna_filtered -a {threads} -m {params.sortmernamem} >> {log} 2>&1
  rm $TMP_R12*
"""

if config["rnafiltering"]["tool"] == "sortmerna":
    rule filter_rna_sortmerna:
        input:
            'Preprocessing/mt.r1.trimmed.fq',
            'Preprocessing/mt.r2.trimmed.fq',
            'Preprocessing/mt.se.trimmed.fq',
            expand(
                "{db}/sortmerna/idx/{files}.{ext}",
                files=config["rnafiltering"]["files"],
                ext=['bursttrie_0.dat', 'kmer_0.dat', 'pos_0.dat', 'stats'],
                db=DBPATH)
        output:
            'Preprocessing/mt.r1.trimmed.rna_filtered.fq',
            'Preprocessing/mt.r2.trimmed.rna_filtered.fq',
            'Preprocessing/mt.se.trimmed.rna_filtered.fq',
            'Preprocessing/mt.r1.trimmed.rna.fq',
            'Preprocessing/mt.r2.trimmed.rna.fq',
            'Preprocessing/mt.se.trimmed.rna.fq'
        threads: 8
        params:
            ref = REF,
            sortmernamem = str(SORTME_MEM),
            runtime = "48:00:00",
            mem = MEMCORE
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/preprocessing_filter_rna.log"
        message:
            "filter_rna: sortmerna: filtering rRNA from mt reads."
        shell:
            FILTER_RNA_SHELL
elif config["rnafiltering"]["tool"] == "bbduk":
    # bbmap repo: https://github.com/BioInfoTools/BBMap/
    # bbduk manual: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
    rule filter_rna_bbduk:
        input:
            r1='Preprocessing/mt.r1.trimmed.fq',
            r2='Preprocessing/mt.r2.trimmed.fq',
            se='Preprocessing/mt.se.trimmed.fq',
            refs=expand(os.path.join(DBPATH, "sortmerna/{files}.fasta"), files=config["rnafiltering"]["files"])
        output:
            # reads unmapped to references
            r1u='Preprocessing/mt.r1.trimmed.rna_filtered.fq',
            r2u='Preprocessing/mt.r2.trimmed.rna_filtered.fq',
            seu='Preprocessing/mt.se.trimmed.rna_filtered.fq',
            # reads mapped to references
            r1m='Preprocessing/mt.r1.trimmed.rna.fq',
            r2m='Preprocessing/mt.r2.trimmed.rna.fq',
            sem='Preprocessing/mt.se.trimmed.rna.fq',
            # stats
            stats_pe='Preprocessing/mt.pe.trimmed.refstats.txt',
            stats_se='Preprocessing/mt.se.trimmed.refstats.txt'
        threads: 8
        params:
            runtime = "12:00:00",
            mem = MEMCORE,
            # from manual: ca. 85% of the machine’s physical memory
            mem2="-Xmx%dg" % round(8 * int(MEMCORE[:-1]) * 0.85),
            # filtering (trade-off between sens and spec, affects memory/runtime)
            filt="k=27 qhdist=1 mkh=1", # (issue #12)
            # incl. rev. complement, quality = auto, do not zip, keep read order
            misc="rcomp=t qin=auto qout=auto ziplevel=1 ordered=t"
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/preprocessing_filter_rna.log"
        message:
            "filter_rna: bbmap/bbduk: filtering rRNA from mt reads."
        shell:
            """
            (
            refs=$(echo \"{input.refs}\" | tr ' ' ',') &&
            bbduk.sh in={input.r1} in2={input.r2} out={output.r1u} out2={output.r2u} outm={output.r1m} outm2={output.r2m} ref=${{refs}} refstats={output.stats_pe} threads={threads} {params.mem2} {params.filt} {params.misc} &&
            bbduk.sh in={input.se} out={output.seu} outm={output.sem} ref=${{refs}} refstats={output.stats_se} threads={threads} {params.mem2} {params.filt} {params.misc}
            ) &> {log}
            """
