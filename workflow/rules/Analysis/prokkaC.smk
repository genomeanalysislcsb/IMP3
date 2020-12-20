rule annotate:
    input:
        'Assembly/%s.assembly.merged.fa' % ASS
    output:
        "Analysis/annotation/annotation.filt.gff",
        "Analysis/annotation/prokka.faa",
        "Analysis/annotation/prokka.fna",
        "Analysis/annotation/prokka.ffn",
        "Analysis/annotation/prokka.fsa",
    threads: getThreads(8)
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    log: "logs/analysis_annotate.log"
    conda: ENVDIR + "/IMP_annotation.yaml"
    message: "annotate: Running prokkaC."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        if [ ! -f $CONDA_PREFIX/db/hmm/HAMAP.hmm.h3m ]; then
          {BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --setupdb
        fi 
	{BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --force --outdir Analysis/annotation --prefix prokka --noanno --cpus {threads} --metagenome {input[0]} >> {log} 2>&1
        
	# Prokka gives a gff file with a long header and with all the contigs at the bottom.  The command below removes the
        # And keeps only the gff table.

        LN=`grep -Hn "^>" Analysis/annotation/prokka.gff | head -n1 | cut -f2 -d ":" || if [[ $? -eq 141 ]]; then true; else exit $?; fi`
        LN1=1
        LN=$(($LN-$LN1))
        head -n $LN Analysis/annotation/prokka.gff | grep -v "^#" | sort | uniq | grep -v "^==" > {output[0]}
        """

