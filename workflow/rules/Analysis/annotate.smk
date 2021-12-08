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
        runtime = "48:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_annotation.yaml"
    threads: getThreads(12)
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
    params:
        config=config["mantis"]
    message:
        "Mantis: config file"
    run:
        with open(output[0], "w+") as ofile:
            # default HMMs
            for ref_name, ref_path in params.config["default"].items():
                ofile.write(f'{ref_name}={ref_path}\n')
            # custom HMMs
            if params.config["custom"]:
                for ref_path in params.config["custom"]:
                    ofile.write(f'custom_ref={ref_path}\n')
            # weights
            for weights_name, weights_value in params.config["weights"].items():
                ofile.write(f'{weights_name}={weights_value}\n')

# Mantis: protein annotation
rule mantis_annot:
    input:
        faa="Analysis/annotation/prokka.faa",
        config="mantis.config",
    output:
        outtab="Analysis/annotation/intermediary/output_annotation.tsv",
        inttab="Analysis/annotation/intermediary/integrated_annotation.tsv",
        contab="Analysis/annotation/intermediary/consensus_annotation.tsv",
        intgff="Analysis/annotation/intermediary/integrated_annotation.gff",
        congff="Analysis/annotation/intermediary/consensus_annotation.gff",
        keggmodules="Analysis/annotation/intermediary/kegg_modules.tsv",
        kosample="Analysis/annotation/intermediary/sample_kos.tsv",
    log:
        "logs/analysis_mantis.log"
    params:
        odir=lambda wildcards, output: os.path.dirname(output[0]),
        params="-gff -km -da heuristic",
    threads: 5
    resources:
        runtime = "24:00:00", # TODO: need runtime estimate for large samples
        mem = MEMCORE,
    message:
        "Mantis: annotate {input}"
    conda:
        os.path.join(ENVDIR, "IMP_mantis.yaml")
    shell:
        "python {MODDIR}/mantis/ run_mantis --target {input.faa} --output_folder {params.odir} --mantis_config {input.config} {params.params} --cores {threads} &> {log}"

# Mantis: reformat consensus output (to be imported in Python/R)
# rule mantis_reformat_consensus:
#     input:
#         "Analysis/annotation/intermediary/mantis/consensus_annotation.tsv"
#     output:
#         temp("Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.tsv")
#     params:
#         runtime = "1:00:00",
#         mem = MEMCORE
#     message:
#         "Mantis: reformat {input}"
#     shell:
#         # merge annotations after "|" and remove "|"
#         "paste <(cut -d '|' -f1 {input} | sed 's/\\t$//') <(cut -d '|' -f2- {input} | sed 's/^\\t//;s/\\t/,/g') > {output}"

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
# rule makegff_mantis:
#     input:
#         gff="Analysis/annotation/annotation.filt.gff",
#         mantis="Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.tsv"
#     output:
#         temp("Analysis/annotation/intermediary/mantis/consensus_annotation.gff")
#     params:
#         runtime = "1:00:00",
#         mem = MEMCORE,
#     threads: 1
#     message:
#         "Mantis: create GFF from {input}"
#     run:
#         import re
#         import numpy
#         import pandas
#         import warnings
#         # read in
#         gff   = pandas.read_csv(input.gff, sep="\t", header=None, names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr"])
#         annot = pandas.read_csv(input.mantis, sep="\t", header=0, index_col=0)
#         annot.replace(to_replace=numpy.nan, value="N/A", inplace=True) # replace None by "N/A" (e.g. no "Links")
#         # to comply with the GFF attributes format:
#         #   replace ";" by "," (";" is used to concat attributes)
#         #   replace "=" by "-" ("=" is used to concat key/value pairs)
#         for cname in ["HMM_Hits", "Links"]:
#             annot[cname] = annot[cname].str.replace(';', ',', regex=True).replace('=', '-', regex=True)
#         # GFF row names = protein ID (extracted from attr column)
#         gff.set_index(keys=gff.attr.map(lambda x: [re.sub("^ID=", "", a) for a in x.split(";") if re.match("ID=", a)][0]), inplace=True)
#         # protein IDs: difference/intersection w/ IDs from Mantis output
#         pids_both = gff.index.intersection(annot.index)
#         pids_diff = gff.index.difference(annot.index)
#         # add Mantis annotations
#         if len(pids_both) > 0:
#             gff.loc[pids_both,"attr"] = gff.loc[pids_both,"attr"].astype(str) + ";mantis_hmmhits=" + annot.loc[pids_both,"HMM_Hits"] + ";mantis_links=" + annot.loc[pids_both,"Links"]
#         else: # raise a warning if no proteins from the GFF have a Mantis annotation
#             warnings.warn("No proteins from %s could be annotated with Mantis (%s" % (input.gff, input.mantis))
#         # add missing Mantis annotations
#         if len(pids_diff) > 0:
#             gff.loc[pids_diff,"attr"] = gff.loc[pids_diff,"attr"].astype(str) + ";mantis_hmmhits=NA;mantis_links=NA"
#         # save
#         gff.to_csv(output[0], sep="\t", header=False, index=False)

# Merge GFFs
rule mergegff:
    input:
        "Analysis/annotation/annotation.filt.gff", # from prokkaC.smk
        expand("Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff", db=config["hmm_DBs"].split()),
        ["Analysis/annotation/intermediary/mantis/consensus_annotation.gff"] if config["annotation"] == "mantis" else []
    output:
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    resources:
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
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    params:
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
