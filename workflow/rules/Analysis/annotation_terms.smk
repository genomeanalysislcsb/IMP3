# One GFF per annotation database from the merged annotation (one line per feature and term; Pfam names -> accessions,
# KEGG pathways once as koNNNNN), for the term counts (scripts/split_annotation_terms.py). The merged annotation holds
# the terms of all sources (Mantis consensus, Bakta), so counts and annotation agree.
PFAM_HMM = config.get("annotation_terms", {}).get("pfam_hmm", "") or os.path.join(DBPATH, "mantis", "pfam", "Pfam-A.hmm")

if "ANNOTATION_TERMS_RULES" not in globals():
    ANNOTATION_TERMS_RULES = True

    rule split_annotation_terms:
        input:
            gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
        output:
            temp("Analysis/annotation/intermediary/terms/annotation.{db}.gff")
        wildcard_constraints:
            # a database name has no dots: otherwise per_contig_counts/<omic>.mantis.<db>.per_contig.counts.tsv
            # could also be made by the total-count rules (db = "<db>.per_contig"), and Snakemake stops as ambiguous
            db = "[A-Za-z0-9_-]+"
        params:
            script = os.path.join(SRCDIR, "split_annotation_terms.py"),
            pfam = ("--pfam-hmm %s" % PFAM_HMM if os.path.exists(PFAM_HMM) else "") + (" --pfam-metadata %s" % os.path.join(os.path.dirname(PFAM_HMM), "metadata.tsv") if os.path.exists(os.path.join(os.path.dirname(PFAM_HMM), "metadata.tsv")) else "")
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            1
        log:
            "logs/analysis/annotation_terms.{db}.log"
        message:
            "split_annotation_terms: {wildcards.db} terms of the merged annotation."
        shell:
            "python {params.script} --gff {input.gff} --db {wildcards.db} {params.pfam} --out {output} > {log} 2>&1"
