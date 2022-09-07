import math

include:
    "mantis_config.smk"

# Mantis: protein annotation
# NOTE: assumming that mantis has been setup
rule mantis_annot:
    input:
        faa="Analysis/annotation/prokka.faa",
        config=rules.mantis_config.output.config_file,
    output:
        outtab="Analysis/annotation/intermediary/mantis/output_annotation.tsv",
        inttab="Analysis/annotation/intermediary/mantis/integrated_annotation.tsv",
        contab="Analysis/annotation/intermediary/mantis/consensus_annotation.tsv",
        intgff="Analysis/annotation/intermediary/mantis/integrated_annotation.gff",
        congff="Analysis/annotation/intermediary/mantis/consensus_annotation.gff",
        keggmodules="Analysis/annotation/intermediary/mantis/kegg_modules.tsv",
        kosample="Analysis/annotation/intermediary/mantis/sample_kos.tsv",
    log:
        "logs/analysis_mantis.log"
    params:
        odir=lambda wildcards, output: os.path.dirname(output.outtab),
        params="-gff -km -da heuristic",
    threads:
        getThreads(128)
    resources:
        runtime="12:00:00",
        mem=MEMCORE
    message:
        "Mantis: annotate {input}"
    conda:
        os.path.join(ENVDIR, "IMP_mantis.yaml")
    shell:
        """
        rm -rf Analysis/annotation/intermediary/mantis*
        mantis run -i {input.faa} \
                   -o {params.odir} \
                   -mc {input.config} \
                   {params.params} \
                   -c {threads} &> {log}
        """

# Mantis: reformat consensus output for easier downstream processing
# NOTE: Cannot use Mantis' GFF file as it does not contain contigs as sequence-records.
#       The GFF file is created from the tabular output.
rule mantis_reformat_consensus:
    input:
        contab=rules.mantis_annot.output.contab
    output:
        contab=temp("Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.tsv")
    resources:
        runtime="1:00:00",
        mem=MEMCORE
    message:
        "Mantis: reformat {input}"
    shell:
        # split by "|", merge annotations (replace "\t" by "|")
        # NOTE: using "|" as it should be very unlikely to have it in the annotation text
        "paste <(cut -d '|' -f1 {input.contab} | sed 's/\\t$//') <(cut -d '|' -f2- {input.contab} | sed 's/^\\t//;s/\\t/|/g') > {output.contab}"

checkpoint makegff_mantis:
    input:
        gff="Analysis/annotation/annotation.filt.gff", # TODO: use rules.annotate.output.XXX when the outputs are named in that rule
        mantis=rules.mantis_reformat_consensus.output.contab,
    output:
        gff=temp("Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.gff")
    params:
        sep="|", # separator used in mantis_reformat_consensus
        prefix="Analysis/annotation/intermediary/mantis/consensus_annotation.mantis"
    resources:
        runtime="12:00:00",
        mem=MEMCORE
    threads:
        # Ensure at least 10gb of mem
        max(int(math.ceil(10 / (int(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]]))), 1)
    message:
        "Mantis: Create GFF from {input}"
    log:
        "logs/analysis_mantis_make_gff.log"
    script:
        os.path.join(SRCDIR, "Analysis", "mantis_gff.py")
