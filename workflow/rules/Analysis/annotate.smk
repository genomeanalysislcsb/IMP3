rule hmmer:
    input:
        "Analysis/annotation/prokka.faa"
    output:
        "Analysis/annotation/intermediary/prokka.faa.{db}.hmmscan"
    params: 
        lambda wildcards: config["hmm_settings"][wildcards.db]["cutoff"],
        dbs = DBPATH + "/hmm/{db}",
        runtime = "48:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_annotation.yaml"
    threads: 12
    log: "logs/analysis_hmmer.{db}.log"
    message: "hmmer: Running HMMER for {wildcards.db}."
    shell:
        """
         hmmsearch --cpu {threads} {params[0]} --noali --notextw \
          --tblout {output} {params.dbs}/*.hmm {input} >/dev/null 2> {log}
        """

rule makegff:
    input:
        "Analysis/annotation/intermediary/prokka.faa.{db}.hmmscan",
        "Analysis/annotation/annotation.filt.gff",
        "Analysis/annotation/prokka.faa"
    output:
        "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    params:
        "{db}",
        lambda wildcards: config["hmm_settings"][wildcards.db]["trim"],
        runtime = "4:00:00",
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_annotation.yaml"
    log: "logs/analysis_makegff.{db}.log"
    message: "makegff: Adding hmmer results of {wildcards.db} to gff."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SRCDIR}/hmmscan_addBest2gff.pl -i {input[0]} -a {input[1]} \
         -n {params[0]} {params[1]} -o {output} -g $(grep ">" {input[2]} | wc -l) > {log} 2>&1
        """

rule mergegff:
    input:
        "Analysis/annotation/annotation.filt.gff",
        expand("Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff",db=config["hmm_DBs"].split())
    output:
        report("Analysis/annotation/annotation_CDS_RNA_hmms.gff",category="Analysis")
    params:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    log: "logs/analysis_mergegff.log"
    message: "mergegff: Merging gffs with hmmer results."
    shell:
        """
        {SRCDIR}/mergegffs.pl {output} {input} > {log} 2>&1
        """

rule tar_annotation:
    input:
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        expand("Analysis/annotation/{type}.{db}.counts.tsv",type=TYPES,db=config["hmm_DBs"].split())
    output:
        "Analysis/annotation/intermediary.tar.gz"
    threads: 1
    params:
        runtime = "8:00:00",
        mem = MEMCORE,
        intermediary = "Analysis/annotation/intermediary/"
    log: "logs/annotation_tar_annotation.log"
    message: "tar_annotation: Compressing intermediary steps of annotation."
    shell:
       """
       tar cvzf {output} {params.intermediary} >> {log} 2>&1 && rm -r {params.intermediary} >> {log} 2>&1
       """

        
        
