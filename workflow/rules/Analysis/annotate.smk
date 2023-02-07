import math

# Rules for annotation w/ Mantis
include:
    "mantis.smk"

# HMMER: protein annotation
rule hmmer:
    input:
        "Analysis/annotation/prokka.faa"
    output:
        "Analysis/annotation/intermediary/prokka.faa.{db}.hmmscan"
    params: 
        lambda wildcards: config["hmm_settings"][wildcards.db]["cutoff"],
        dbs = DBPATH + "/hmm/{db}"
    resources:
        runtime = "12:00:00",
        mem_mb = MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    threads: getThreads(12)
    log: "logs/analysis_hmmer.{db}.log"
    message: "hmmer: Running HMMER for {wildcards.db}."
    shell:
        """
         hmmsearch --cpu {threads} \
                   {params[0]} \
                   --noali \
                   --notextw \
                   --tblout {output} \
                   {params.dbs}/*.hmm \
                   {input} >/dev/null 2> {log}
        """

# GFF from HMMER output
rule makegff:
    input:
        "Analysis/annotation/intermediary/prokka.faa.{db}.hmmscan",
        "Analysis/annotation/annotation.filt.gff",
        "Analysis/annotation/prokka.faa"
    output:
        "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    params:
        "{db}",
        lambda wildcards: config["hmm_settings"][wildcards.db]["trim"]
    resources:
        runtime = "02:00:00",
        mem_mb = MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/analysis_makegff.{db}.log"
    message:
        "makegff: Adding hmmer results of {wildcards.db} to gff."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SRCDIR}/hmmscan_addBest2gff.pl -i {input[0]} \
                                        -a {input[1]} \
                                        -n {params[0]} {params[1]} \
                                        -o {output} \
                                        -g $(grep ">" {input[2]} | wc -l) > {log} 2>&1
        """

# Merge GFFs
rule mergegff:
    input:
        "Analysis/annotation/annotation.filt.gff", # from prokkaC.smk
        expand("Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff", db=config["hmm_DBs"].split())
          if config["annotation"] == "hmmer" else [],
        # ["Analysis/annotation/intermediary/annotation.CDS.RNA.essential.gff"] if (config["annotation"] != "hmmer" and "binny" in BIN_STEPS) else [],
        [],  ### ???
        [rules.makegff_mantis.output] if config["annotation"] == "mantis" else []
    output:
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    resources:
        runtime = "04:00:00",
        mem_mb = MEMCORE_MB
    threads:
        # Ensure at least 10gb of mem
        max(int(math.ceil(10 / (int(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]]))), 1)
    log:
        "logs/analysis_mergegff.log"
    message:
        "mergegff: Merging gffs with hmmer results."
    shell:
        """
        {SRCDIR}/mergegffs.pl {output} {input} > {log} 2>&1
        """

# Compress intermediate files
rule tar_annotation:
    input:
        # Make sure its done at the very end
        "status/workflow.done",
        gff="Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        # from featureCount.smk
        counts_hmmer=expand("Analysis/annotation/{type}.hmmer.{db}.counts.tsv", type=TYPES,
                            db=config["hmm_DBs"].split()) if config["annotation"] == "hmmer" else [],
        counts_mantis=aggregate_mantis_fc if config["annotation"] == "mantis" else [],
        bins="Binning/MetaWrap/bins.stats" if 'binning' in IMP_STEPS else []
    output:
        "Analysis/annotation/intermediary.tar.gz"
    threads:
        # Ensure at least 10gb of mem
        max(int(math.ceil(10 / (int(MEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]]))), 8)
    resources:
        runtime = "04:00:00",
        mem_mb = MEMCORE_MB
    params:
        intermediary = "Analysis/annotation/intermediary/"
    log:
        "logs/annotation_tar_annotation.log"
    message:
        "tar_annotation: Compressing intermediary steps of annotation."
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    shell:
       """
       tar -cv {params.intermediary} | pigz --stdout \
                                            --best \
                                            --processes {threads} \
                                            --rsyncable \
                                            --verbose > {output} 2> {log} && \
       rm -rf {params.intermediary} >> {log} 2>&1
       """
