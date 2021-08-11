rule cut_rRNA:
    input:
        "Assembly/%s.assembly.merged.fa" %ASS,
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "Binning/binny/%s.assembly.merged.cut.fa" %ASS
    log: "logs/binning_cut_rRNA.log"
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_annotation.yaml" 
    message: "cut_rRNA: Cutting contigs for vizbin."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SRCDIR}/fastaExtractCutRibosomal1000.pl -f {input[0]} -g {input[1]} -l {log} -o {output} -c {config[binning][vizbin][cutoff]}
        """

rule vizbin:
    input:
        "Binning/binny/%s.assembly.merged.cut.fa" %ASS
    output:
        "Binning/binny/%s.vizbin.with-contig-names.points" %ASS
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    params:
        Xmx = str(round(BIGMEMTOTAL*0.75)) + "g"
    conda: ENVDIR + "/IMP_annotation.yaml" 
    threads: getThreads(BIGCORENO) 
    log: "logs/binning_vizbin.%s.log" %ASS
    message: "vizbin: Running VizBin."
    shell:
        """
        TMP_VIZBIN=$(mktemp --tmpdir=Analysis -dt "VIZBIN_XXXXXX")
        java -Xmx {params.Xmx} -jar {BINDIR}/VizBin-dist.jar \
         -a {config[binning][vizbin][dimension]} \
         -c {config[binning][vizbin][cutoff]} \
         -i {input} \
         -o $TMP_VIZBIN/data.points \
         -k {config[binning][vizbin][kmer]} \
         -p {config[binning][vizbin][perp]} > {log} 2>&1

         if [ -f $TMP_VIZBIN/data.points ]
           then
             paste <(grep "^>" {input} | sed -e 's/>//') \
              <(cat $TMP_VIZBIN/data.points | sed -e 's/,/\t/') > {output}
           else
             touch {output} status/binning_vizbin.impossible
           fi
        rm -rf $TMP_VIZBIN
        """
        
rule vizbin_mt:
    input:
        "Assembly/mt.assembly.merged.fa"
    output:
        "Binning/binny/mt.vizbin.with-contig-names.points"
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    params:
        Xmx=str(round(BIGMEMTOTAL*0.75)) + "g"
    threads: getThreads(BIGCORENO)
    log: "logs/binning_vizbin.mt.log"
    message: "vizbin: Running VizBin."
    shell:
        """
        TMP_VIZBIN=$(mktemp --tmpdir=Analysis -dt "VIZBIN_XXXXXX")
        java -Xmx{params.Xmx} -jar {BINDIR}/VizBin-dist.jar \
         -a {config[binning][vizbin][dimension]} \
         -c {config[binning][vizbin][cutoff]} \
         -i {input} \
         -o $TMP_VIZBIN/data.points \
         -k {config[binning][vizbin][kmer]} \
         -p {config[binning][vizbin][perp]} > {log} 2>&1

         if [ -f $TMP_VIZBIN/data.points ]
           then
             paste <(grep "^>" {input} | sed -e 's/>//') \
              <(cat $TMP_VIZBIN/data.points | sed -e 's/,/\t/') > {output}
           else
             touch {output} status/binning_vizbin.impossible
           fi
        rm -rf $TMP_VIZBIN
        """
