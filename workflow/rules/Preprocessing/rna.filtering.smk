# https://github.com/biocore/sortmerna/
# https://github.com/biocore/sortmerna/wiki/User-manual-v4.0
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
        # list of reference FASTA files
        ref = ' '.join([
            '--ref %s' % (a) for a in expand(
                "{path}/sortmerna/{files}.fasta",
                files=config["sortmerna"]["files"],
                path=DBPATH
            )
        ]),
        # --fastx: output aligned reads
        # -m: memory requirements for indexing
        params = "--fastx -m {mem}".format(mem=int(MEMCORE[:-1]) * 1000),
        # --otu-map: OTU map (input to QIIME's make_otu_table.py)
        # --paired_in: flags the PE reads as "aligned", when either of them is aligned
        params_pe = "--otu_map --paired_in",
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log:
        "logs/preprocessing_filter_rna.log"
    message:
        "filter_rna: filtering rRNA from mt reads."
    shell:
        # NOTE:
        # sortmerna adds an extension to output files automatically
        # sub-folder kvdb/ in the work directory has to be cleaned before running sortmerna
        """
        IDX_DIR=$(dirname $(dirname {input.idx[0]}))
        TMP_R12=$(mktemp --tmpdir={TMPDIR} -t "tmp_XXXXXX.R12")
        FQ_SE_FILT={output.filt_se}
        FQ_SE_RRNA={output.rrna_se}
        if [ -d {TMPDIR}/kvdb ]; then
            rm -Rf {TMPDIR}/kvdb
        fi
        
        sortmerna {params.ref} --idx-dir ${{IDX_DIR}} --reads {input.r1} --reads {input.r2} \
        --workdir {TMPDIR} --aligned {TMPDIR}/${{TMP_R12}}.rrna --other {TMPDIR}/${{TMP_R12}}.rrna_filt \
        --threads {threads} {params.params} {params.params_pe} &> {log}
        
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
        sortmerna {params.ref} --idx-dir ${{IDX_DIR}} --reads {input.se} \
        --workdir {TMPDIR} --aligned ${{FQ_SE_RRNA%.*}} --other ${{FQ_SE_FILT%.*}} \
        --threads {threads} {params.params} &>> {log}
        rm -Rf {TMPDIR}/${{TMP_R12}}* 
        """
