FASTAFILES = expand("{path}/{files}.fasta", files=config["preprocessing"]["filtering"]["rrna"], path=DBPATH + "/sortmerna")

REF = ' '.join('--ref %s' % (a) for a in zip(FASTAFILES))

SORTME_MEM = get_mem_norm() * 1000

FILTER_RNA_SHELL = """
  TMP_R12=$(mktemp --tmpdir=${{TMPDIR}} -t "tmp_XXXXXX.R12")
  if [ -d ${{TMPDIR}}/kvdb ]; then rm -Rf ${{TMPDIR}}/kvdb; fi 
  sortmerna {params.ref} --reads {input[0]} --reads {input[1]} --workdir tmp --fastx --aligned ${{TMPDIR}}/$TMP_R12.rrna \
   --other ${{TMPDIR}}/$TMP_R12.rrna_filt --threads {threads} --idx-dir {DBPATH}/sortmerna --otu_map --paired_in -m {params.sortmernamem} >> {log} 2>&1
  otufile=`find . -name otu_map.txt`
  touch Preprocessing/otu_map.txt
  cp $otufile Preprocessing/
  if [ -s $otufile ]; then echo "# WARNING: OTU from sortmerna file does not exists" >status/sortmerna.warning; fi
  {SRCDIR}/unmerge-paired-reads.sh ${{TMPDIR}}/$TMP_R12.rrna.fq {output[3]} {output[4]} >> {log} 2>&1
  {SRCDIR}/unmerge-paired-reads.sh ${{TMPDIR}}/$TMP_R12.rrna_filt.fq {output[0]} {output[1]} >> {log} 2>&1
  rm -f ${{TMPDIR}}/kvdb/*
  sortmerna {params.ref} --reads {input[2]} --fastx --aligned Preprocessing/mt.se.trimmed.rna --idx-dir {DBPATH}/sortmerna --workdir ${{TMPDIR}} \
   --other Preprocessing/mt.se.trimmed.rna_filtered --threads {threads} -m {params.sortmernamem} >> {log} 2>&1
  rm -Rf ${{TMPDIR}}/$TMP_R12* 
"""

rule filter_rna:
    input:
        'Preprocessing/mt.r1.trimmed.fq',
        'Preprocessing/mt.r2.trimmed.fq',
        'Preprocessing/mt.se.trimmed.fq',
        expand(
            "{db}/sortmerna/idx/{files}.{ext}",
            files=config["preprocessing"]["filtering"]["rrna"],
            ext=['bursttrie_0.dat', 'kmer_0.dat', 'pos_0.dat', 'stats'],
            db=DBPATH)
    output:
        'Preprocessing/mt.r1.trimmed.rna_filtered.fq',
        'Preprocessing/mt.r2.trimmed.rna_filtered.fq',
        'Preprocessing/mt.se.trimmed.rna_filtered.fq',
        'Preprocessing/mt.r1.trimmed.rna.fq',
        'Preprocessing/mt.r2.trimmed.rna.fq',
        'Preprocessing/mt.se.trimmed.rna.fq',
	'Preprocessing/otu_map.txt'
    threads: 8
    resources:
        runtime = "48:00:00",
        mem = get_mem_norm() * 8
    params:
        ref = REF,
        sortmernamem = str(SORTME_MEM)
    conda: ENVDIR + "/IMP_preprocessing.yaml"
    log: "logs/preprocessing_filter_rna.log"
    message: "filter_rna: filtering rRNA from mt reads."
    shell:
        FILTER_RNA_SHELL
