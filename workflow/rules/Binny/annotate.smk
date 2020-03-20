rule hmmer_essential:
    input:
        "Analysis/annotation/prokka.faa"
    output:
        "Analysis/annotation/intermediary/prokka.faa.essential.hmmscan"
    params: 
        dbs = DBPATH + "/hmm/essential",
        runtime = "48:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_annotation.yaml"
    threads: 12
    log: "logs/analysis_hmmer.essential.log"
    message: "hmmer: Running HMMER for essential."
    shell:
        """
         hmmsearch --cpu {threads} --cut_tc --noali --notextw \
          --tblout {output} {params.dbs}/*.hmm {input} >/dev/null 2> {log}
        """

rule makegff:
    input:
        "Analysis/annotation/intermediary/prokka.faa.essential.hmmscan",
        "Analysis/annotation/annotation.filt.gff",
        "Analysis/annotation/prokka.faa"
    output:
        "Analysis/annotation/intermediary/annotation.CDS.RNA.essential.gff"
    params:
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
         -n essential -o {output} -g $(grep ">" {input[2]} | wc -l) > {log} 2>&1
        """

rule mergegff:
    input:
        "Analysis/annotation/annotation.filt.gff",
        "Analysis/annotation/intermediary/annotation.CDS.RNA.essential.gff"
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
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
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
       tar cvzf {output} {params.intermediary} && rm -r {params.intermediary} >> {log} 2>&1
       """

        
        
