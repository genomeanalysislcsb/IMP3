rule filter_rna:
    input:
        r1='Preprocessing/mt.r1.trimmed.fq',
        r2='Preprocessing/mt.r2.trimmed.fq',
        se='Preprocessing/mt.se.trimmed.fq',
        idx=expand(
            "{db}/sortmerna/idx/{files}.{ext}",
            db=DBPATH,
            files=config["sortmerna"]["files"],
            ext=['bursttrie_0.dat', 'kmer_0.dat', 'pos_0.dat', 'stats']
        )
    output:
        filt_r1='Preprocessing/mt.r1.trimmed.rna_filtered.fq',
        filt_r2='Preprocessing/mt.r2.trimmed.rna_filtered.fq',
        filt_se='Preprocessing/mt.se.trimmed.rna_filtered.fq',
        rrna_r1='Preprocessing/mt.r1.trimmed.rna.fq',
        rrna_r2='Preprocessing/mt.r2.trimmed.rna.fq',
        rrna_se='Preprocessing/mt.se.trimmed.rna.fq',
	    otu='Preprocessing/otu_map.txt'
    threads: getThreads(8)
    resources:
        runtime = "48:00:00",
        mem = MEMCORE
    params:
        ref = ' '.join([
            '--ref %s' % (a) for a in expand(
                "{path}/sortmerna/{files}.fasta",
                files=config["sortmerna"]["files"],
                path=DBPATH
            )
        ]),
        sortmernamem = str(int(MEMCORE[:-1]) * 1000)
    conda: ENVDIR + "/IMP_preprocessing.yaml"
    log: "logs/preprocessing_filter_rna.log"
    message: "filter_rna: filtering rRNA from mt reads."
    shell:
        """
        IDX_DIR=$(dirname $(dirname {input.idx[0]}))
        TMP_R12=$(mktemp --tmpdir={TMPDIR} -t "tmp_XXXXXX.R12")
        FQ_SE_FILT={output.filt_se}
        FQ_SE_RRNA={output.rrna_se}
        if [ -d {TMPDIR}/kvdb ]; then
            rm -Rf {TMPDIR}/kvdb
        fi
        
        sortmerna {params.ref} --reads {input.r1} --reads {input.r2} --workdir {TMPDIR} --fastx --aligned {TMPDIR}/${{TMP_R12}}.rrna \
        --other {TMPDIR}/${{TMP_R12}}.rrna_filt --threads {threads} --idx-dir ${{IDX_DIR}} --otu_map --paired_in -m {params.sortmernamem} &> {log}
        otufile=`find . -name otu_map.txt`

        if [ -s ${{otufile}} ]; then
            cp ${{otufile}} {output.otu}
        else
            touch {output.otu}
            echo "# WARNING: OTU from sortmerna file does not exists" > status/sortmerna.warning
        fi

        {BINDIR}/unmerge-paired-reads.sh {TMPDIR}/${{TMP_R12}}.rrna.fq {output.rrna_r1} {output.rrna_r2} &>> {log}
        {BINDIR}/unmerge-paired-reads.sh {TMPDIR}/${{TMP_R12}}.rrna_filt.fq {output.filt_r1} {output.filt_r2} &>> {log}

        rm -f {TMPDIR}/kvdb/*
        sortmerna {params.ref} --reads {input.se} --fastx --aligned ${{FQ_SE_RRNA%.*}} --idx-dir ${{IDX_DIR}} --workdir {TMPDIR} \
        --other ${{FQ_SE_FILT%.*}} --threads {threads} -m {params.sortmernamem} &>> {log}
        rm -Rf {TMPDIR}/${{TMP_R12}}* 
        """
