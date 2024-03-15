# include rules for the Taxonomy part of the workflow

include:
    "../../Taxonomy/mOTUlink.smk"
include:
    "../../Taxonomy/mOTUs.smk"
include:
    "../../Taxonomy/kraken_reads.smk"
include:
    "../../Taxonomy/kraken_reads_SE.smk"
include:
    "../../Taxonomy/kraken_contigs.smk"
include:
    "../../Taxonomy/gtdbtk_bins.smk"
if config['eukdetect']['run_eukdetect']:
    include:
        "../../Taxonomy/eukdetect.smk"

# master command
if MG or MT:
    # Since we now combine both reports we need to track only the report here
    tax_ctrl = expand("Analysis/taxonomy/kraken/{type}.kraken.{outfile}",type=TYPES,outfile="report")
    tax_ctrl.append(expand("status/{type}.bracken.done",type=TYPES)),
    tax_ctrl.append(expand("Analysis/taxonomy/mOTUs/{type}.mOTU.{out_type}.tsv", type=TYPES, out_type=["counts", "rel_abun"]))
    if config['eukdetect']['run_eukdetect']:
        tax_ctrl.append("status/eukdetect.done")
if "assembly" in IMP_STEPS:
    tax_ctrl.append("Analysis/taxonomy/kraken/%s.kraken.contig.output.gz" % ASS)
    if "analysis" in IMP_STEPS:
        tax_ctrl.append("Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv")
if CONTIGS:
    tax_ctrl = ["Analysis/taxonomy/kraken/%s.kraken.contig.output.gz" % ASS]
    if MG or MT:
        tax_ctrl.append(expand("Analysis/taxonomy/kraken/{type}.kraken.{outfile}",type=TYPES,outfile=["output.gz","report"]))
        if "analysis" in IMP_STEPS:
            tax_ctrl.append("Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv")
if "binning" in IMP_STEPS:
    tax_ctrl.append("status/taxonomy_GTDBtk.done")

rule TAXONOMY:
    input:
        tax_ctrl
    output:
        touch('status/taxonomy.done')
