#
# Taxonomy.smk
#
# This file includes all rules for the taxonomy part of the workflow.
# It contains logic to switch between the standard single-database Kraken
# workflow and the experimental dual-database (GTDB + NT) workflow.
#

# --- Workflow Switching Logic ---
# Check the config for the 'kraken_dual_db' section.
# The 'get' method with a default empty dict prevents errors if the section is missing.
use_dual_db = config.get("kraken_dual_db", {}).get("run_dual_db", False)

# --- Include appropriate rule files based on the selected workflow ---
if use_dual_db:
    # --- EXPERIMENTAL DUAL-DATABASE WORKFLOW ---
    # This path is executed if 'run_dual_db' is true in the config.
    print("INFO: Using experimental dual-database Kraken workflow.")
    # First, include the unified base rules for kraken and processing steps.
    # This file must be included before the rules that inherit from it.
    include: "../../Taxonomy/kraken_base.smk"
    # Then, include the specific dual-db workflow rules that use inheritance.
    include: "../../Taxonomy/kraken_reads_dual.smk"
    include: "../../Taxonomy/kraken_contigs_dual.smk"
else:
    # --- STANDARD SINGLE-DATABASE WORKFLOW ---
    # This is the default path, which uses the original, self-contained rules.
    print("INFO: Using standard single-database Kraken workflow.")
    if len(MG) == 2 or len(MT) == 2:
        include: "../../Taxonomy/kraken_reads.smk"
    else:
        include: "../../Taxonomy/kraken_reads_SE.smk"
    include: "../../Taxonomy/kraken_contigs.smk"

# --- Include other standard taxonomy rules ---
include: "../../Taxonomy/mOTUlink.smk"
include: "../../Taxonomy/mOTUs.smk"
include: "../../Taxonomy/gtdbtk_bins.smk"

# --- Master Rule: TAXONOMY ---
# This rule aggregates all the final outputs from the taxonomy module.
# The 'tax_ctrl' list is built dynamically based on the running workflow.
# Since the final output files are named identically in both workflows,
# this logic does not need to change.

tax_ctrl = []

if MG or MT:
    tax_ctrl.extend(expand("Analysis/taxonomy/kraken/{type}.kraken.report", type=TYPES))
    tax_ctrl.append(expand("status/{type}.bracken.done", type=TYPES))
    tax_ctrl.append(expand("Analysis/taxonomy/mOTUs/{type}.mOTU.{out_type}.tsv", type=TYPES, out_type=["counts", "rel_abun"]))

if "assembly" in IMP_STEPS:
    tax_ctrl.append("Analysis/taxonomy/kraken/%s.kraken.contig.report" % ASS)
    if "analysis" in IMP_STEPS:
        tax_ctrl.append("Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv")

if CONTIGS:
    tax_ctrl = ["Analysis/taxonomy/kraken/%s.kraken.contig.report" % ASS]
    if MG or MT:
        tax_ctrl.append(expand("Analysis/taxonomy/kraken/{type}.kraken.report", type=TYPES))
        if "analysis" in IMP_STEPS:
            tax_ctrl.append("Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv")

if "binning" in IMP_STEPS:
    tax_ctrl.append("status/taxonomy_GTDBtk.done")

rule TAXONOMY:
    input:
        tax_ctrl
    output:
        touch('status/taxonomy.done')
    message:
        "Taxonomy module completed."
