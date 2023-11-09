rule antismash_filter_ass:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS
    output:
        "Assembly/%s.assembly.merged.length_filtered.fa" % ASS
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_annotation.yaml"
    threads: 1
    log: "logs/antismash_filter_ass.log"
    message: "antismash: Filtering assembly."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        cat {input} | perl {SRCDIR}/fasta_filter_length.pl 10000 > {output} 2> {log}
        """

rule antismash_filter_gff:
    input:
        in_gff="Analysis/annotation/annotation.filt.gff",
        in_ass="Assembly/%s.assembly.merged.length_filtered.fa" % ASS
    output:
        "Analysis/annotation/annotation.filt.ass_length.gff"
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    conda: os.path.join(ENVDIR, "IMP_antismash.yml")
    threads: 1
    log: "logs/antismash_filter_gff.log"
    message: "antismash: Filtering gff."
    script:
        os.path.join(SRCDIR, "antismash_gff_filtering.py")

rule antismash:
    input:
        gff="Analysis/annotation/annotation.filt.ass_length.gff",
        assembly="Assembly/%s.assembly.merged.length_filtered.fa" % ASS
    output:
        "Analysis/antiSMASH/index.html"
    params:
        out_dir="Analysis/antiSMASH",
        db=os.path.join(DBPATH,"antismash/databases")
    resources:
        runtime = "4:00:00",
        mem = MEMCORE
    threads: getThreads(8)
    conda: os.path.join(ENVDIR, "IMP_antismash.yml")
    log: "logs/antismash_main.log"
    message: "antismash: Running antiSMASH."
    shell:
        """
        if [ ! -f {DBPATH}/antismash/databases/clusterblast/proteins.dmnd ]; then
           download-antismash-databases --database-dir {DBPATH}/antismash/databases
        fi
        antismash {input.assembly} --output-dir {params.out_dir} --genefinding-gff3 {input.gff} -c {threads} --allow-long-headers --asf --clusterhmmer --cb-knownclusters --rre --cb-subcluster
        """


rule gbk_parsing:
    input:
        "Analysis/antiSMASH/index.html"
    output:
        "Analysis/antiSMASH/antismash.all.gff",
        "Analysis/antiSMASH/antismash.db_xref.gff",
        "Analysis/antiSMASH/antismash.gene_functions.gff",
        "Analysis/antiSMASH/antismash.aSDomain.gff",
        "Analysis/antiSMASH/antismash.MiBig.gff"
    params:
        aS_dir="Analysis/antiSMASH" 
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    conda: os.path.join(ENVDIR, "IMP_antismash.yml")
    log: "logs/antismash_gbk2gff.log"
    message: "antismash: convert results from gbk to gff."
#    script:
#        os.path.join(SRCDIR,"antismash_gbk_parsing.py")
    shell:
        """
        if ls {params.aS_dir}/*contig*gbk 1> /dev/null 2>&1; then
            echo "parsing .gbk output" &> {log}
            python3 {SRCDIR}/antismash_gbk_parsing_shell.py {params.aS_dir} &>> {log}
        else
            echo "no antiSMASH output" &>> {log}
            touch {output}
        fi
        """

        
        
