import math

# Rules for annotation w/ Mantis
include:
    "mantis.smk"

# Merge GFFs
rule mergegff:
    input:
        "Analysis/annotation/bakta/annotation.filt.gff", # from bakta
        rules.makegff_mantis.output
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
        "mergegff: Merging the Bakta and Mantis gffs."
    shadow: "copy-minimal"
    shell:
        """
        {SRCDIR}/mergegffs.pl {output} {input} > {log} 2>&1
        """
