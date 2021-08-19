# include rules for the Taxonomy part of the workflow

include:
    "../../Taxonomy/mOTUlink.smk"
include:
    "../../Taxonomy/mOTUs.smk"
include:
    "../../Taxonomy/kraken_reads.smk"
include:
    "../../Taxonomy/kraken_contigs.smk"
include:
    "../../Taxonomy/gtdbtk_bins.smk"
if config['eukdetect']['run_eukdetect']:
    include:
        "../../Taxonomy/eukdetect.smk"

# master command
if MG or MT:
    tax_ctrl = expand(
        "Analysis/taxonomy/kraken/{type}.{kname}.{outfile}",
        type=TYPES, kname=['kraken', 'se.kraken'], outfile=["output.gz", "report"]
    )
    tax_ctrl.append(expand("Analysis/taxonomy/mOTUs/{type}.mOTU.counts.tsv",type=TYPES))
    if config['eukdetect']['run_eukdetect']:
        tax_ctrl.append("status/eukdetect.done")
if "assembly" in IMP_STEPS:
    tax_ctrl.append("Analysis/taxonomy/kraken/%s.kraken.parsed.tsv" % ASS)
    if "analysis" in IMP_STEPS:
        tax_ctrl.append("Analysis/taxonomy/mOTU_links/intermediary.tar.gz")
if CONTIGS:
    tax_ctrl = ["Analysis/taxonomy/kraken/%s.kraken.parsed.tsv" % ASS]
    if MG or MT:
        tax_ctrl.append(expand(
            "Analysis/taxonomy/kraken/{type}.{kname}.{outfile}",
            type=TYPES, kname=['kraken', 'se.kraken'], outfile=["output.gz", "report"]
        ))
        if "analysis" in IMP_STEPS:
            tax_ctrl.append("Analysis/taxonomy/mOTU_links/intermediary.tar.gz")
if "binning" in IMP_STEPS:
    tax_ctrl.append("status/taxonomy_GTDBtk.done")

rule TAXONOMY:
    input:
        tax_ctrl
    output:
        touch('status/taxonomy.done')
