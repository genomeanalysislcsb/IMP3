import math

# Rules for annotation w/ Mantis
include:
    "mantis.smk"

# HMMER: protein annotation
rule hmmer:
    input:
        "Analysis/annotation/bakta/annotation.faa"
    output:
        "Analysis/annotation/intermediary/annotation.faa.{db}.hmmscan"
    params: 
        lambda wildcards: config["hmm_settings"][wildcards.db]["cutoff"],
        dbs = DBPATH + "/hmm/{db}"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    threads: getThreads(12)
    log: "logs/analysis/hmmer.{db}.log"
    message: "hmmer: Running HMMER for {wildcards.db}."
    shadow: "copy-minimal"
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
        "Analysis/annotation/intermediary/annotation.faa.{db}.hmmscan",
        "Analysis/annotation/bakta/annotation.filt.gff",
        "Analysis/annotation/bakta/annotation.faa"
    output:
        "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    params:
        "{db}",
        lambda wildcards: config["hmm_settings"][wildcards.db]["trim"]
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/analysis/makegff.{db}.log"
    message:
        "makegff: Adding hmmer results of {wildcards.db} to gff."
    shadow: "copy-minimal"
    shell:
        """
        export PATH=${{CONDA_PREFIX}}/bin:$PATH
        export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl
        {SRCDIR}/hmmscan_addBest2gff.pl -i {input[0]} \
                                        -a {input[1]} \
                                        -n {params[0]} {params[1]} \
                                        -o {output} \
                                        -g $(grep ">" {input[2]} | wc -l) > {log} 2>&1
        """

# Merge GFFs
rule mergegff:
    input:
        "Analysis/annotation/bakta/annotation.filt.gff", # from bakta
        expand("Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff", db=config["hmm_DBs"].split())
          if config["annotation"] == "hmmer" else [],
        # ["Analysis/annotation/intermediary/annotation.CDS.RNA.essential.gff"] if (config["annotation"] != "hmmer" and "binny" in BIN_STEPS) else [],
        [],  ### ???
        [rules.makegff_mantis.output] if config["annotation"] == "mantis" else []
    output:
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        # Ensure at least 50gb of mem
        max(int(math.ceil(50 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    log:
        "logs/analysis/mergegff.log"
    message:
        "mergegff: Merging gffs with hmmer results."
    shadow: "copy-minimal"
    shell:
        """
        {SRCDIR}/mergegffs.pl {output} {input} > {log} 2>&1
        """
