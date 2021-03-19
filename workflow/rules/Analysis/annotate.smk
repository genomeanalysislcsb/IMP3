# HMMER: protein annotation
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

# Mantis: create config file from IMP config
localrules: mantis_config
rule mantis_config:
    output:
        "mantis.config"
    message:
        "Mantis: create config"
    run:
        with open(output[0], "w") as ofile:
            # default HMMs
            for hmm_name, hmm_path in config["mantis"]["default"].items():
                ofile.write("%s=%s\n" % (hmm_name, hmm_path))
            # custom HMMs
            for hmm_path in config["mantis"]["custom"]:
                ofile.write("custom_hmm=%s\n" % hmm_path)
            # weights
            for weights_name, weights_value in config["mantis"]["weights"].items():
                ofile.write("%s=%f\n" % (weights_name, weights_value))

# Mantis: protein annotation
# NOTE: check installation before use: python submodules/mantis/ check_installation
rule mantis:
    input:
        faa="Analysis/annotation/prokka.faa",
        config="mantis.config"
    output:
        outtab="Analysis/annotation/intermediary/mantis/output_annotation.tsv",
        inttab="Analysis/annotation/intermediary/mantis/integrated_annotation.tsv",
        contab="Analysis/annotation/intermediary/mantis/consensus_annotation.tsv"
    log:
        "logs/analysis_mantis.log"
    params:
        runtime = "24:00:00", # TODO: need runtime estimate for large samples (test: < 2h w/ 10 cores)
        mem = MEMCORE,
        # --default_workers, --chunk_size: use default
        # --hmmer_threads 1: using more threads does not really improve the runtime and consumes more RAM
        params="--hmmer_threads 1",
        odir=lambda wildcards, output: os.path.dirname(output[0])
    threads: 12
    conda:
        os.path.join(ENVDIR, "IMP_MANTIS.yaml")
    message:
        "Mantis: annotate {input}"
    shell:
        "(if [ -d {params.odir} ]; then ls -l {params.odir}; fi && python {MODDIR}/mantis/ run_mantis --target {input.faa} --output_folder {params.odir} --mantis_config {input.config} {params.params}) &> {log}"

# Mantis: reformat consensus output (to be imported in Python/R)
rule mantis_reformat_consensus:
    input:
        "Analysis/annotation/intermediary/mantis/consensus_annotation.tsv"
    output:
        temp("Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.tsv")
    params:
        runtime = "1:00:00",
        mem = MEMCORE
    message:
        "Mantis: reformat {input}"
    shell:
        # merge annotations after "|", remove "|" separator
        "paste <(cut -d '|' -f1 {input} | sed 's/\\t$//') <(cut -d '|' -f2 {input} | sed 's/^\\t//;s/\\t/;/g') > {output}"

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

# Mantis: GFF from (reformatted) consensus output
rule makegff_mantis:
    input:
        gff="Analysis/annotation/annotation.filt.gff",
        mantis="Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.tsv"
    output:
        temp("Analysis/annotation/intermediary/mantis/consensus_annotation.gff")
    params:
        runtime = "1:00:00",
        mem = MEMCORE,
    threads: 1
    message:
        "Mantis: create GFF from {input}"
    run:
        import re
        import pandas
        # read in
        gff   = pandas.read_csv(input.gff, sep="\t", header=None, names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr"])
        annot = pandas.read_csv(input.mantis, sep="\t", header=0, index_col=0)
        # keep only entries w/ annotation
        gff   = gff.loc[gff.index.intersection(annot.index),:]
        # GFF row names: extracted from attr column
        gff.set_index(keys=gff.attr.map(lambda x: [re.sub("^ID=", "", a) for a in x.split(";") if re.match("ID=", a)][0]), inplace=True)
        # add annotations
        gff.loc[:,"attr"] = gff.loc[:,"attr"].astype(str) + ";" + annot.loc[gff.index,"Links"]
        # save
        gff.to_csv(output[0], sep="\t", header=False, index=False)

# Merge GFFs
rule mergegff:
    input:
        "Analysis/annotation/annotation.filt.gff", # from prokkaC.smk
        expand("Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff", db=config["hmm_DBs"].split()),
        ["Analysis/annotation/intermediary/mantis/consensus_annotation.gff"] if config["annotation"] == "mantis" else []
    output:
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
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

# Compress intermediate files
rule tar_annotation:
    input:
        gff="Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        # from featureCount.smk
        counts=expand("Analysis/annotation/{type}.{db}.counts.tsv", type=TYPES, db=config["hmm_DBs"].split())
    output:
        "Analysis/annotation/intermediary.tar.gz"
    threads: 1
    params:
        runtime = "8:00:00",
        mem = MEMCORE,
        intermediary = "Analysis/annotation/intermediary/"
    log:
        "logs/annotation_tar_annotation.log"
    message:
        "tar_annotation: Compressing intermediary steps of annotation."
    shell:
       """
       tar cvzf {output} {params.intermediary} >> {log} 2>&1 && \
       rm -r {params.intermediary} >> {log} 2>&1
       """