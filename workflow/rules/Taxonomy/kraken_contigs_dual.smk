#
# kraken_contigs_dual.smk
#
# This file contains the complete Snakemake rules for the experimental
# dual-database Kraken2 workflow for assembled contigs, using rule inheritance.
#

# --- Define Kraken2 runs for contigs against each database ---

use rule kraken2_base as kraken_contigs_gtdb with:
    input: "Assembly/%s.assembly.merged.fa" % ASS
    output: temp("Analysis/taxonomy/kraken/%s.gtdb.kraken.contig.output.gz" % ASS)
    params: database = os.path.join(DBPATH, config["kraken_dual_db"]["gtdb_name"]), extra_opts = ""
    log: "logs/analysis/taxonomy/kraken_contigs_dual.gtdb.log"
    message: "Dual-DB Kraken (GTDB): Running kraken2 on contigs."

use rule kraken2_base as kraken_contigs_nt with:
    input: "Assembly/%s.assembly.merged.fa" % ASS
    output: temp("Analysis/taxonomy/kraken/%s.nt.kraken.contig.output.gz" % ASS)
    params: database = os.path.join(DBPATH, config["kraken_dual_db"]["nt_name"]), extra_opts = ""
    log: "logs/analysis/taxonomy/kraken_contigs_dual.nt.log"
    message: "Dual-DB Kraken (NT): Running kraken2 on contigs."

# --- Combine results and generate final report ---

use rule combine_kraken_base as combine_kraken_contigs with:
    input: gtdb="Analysis/taxonomy/kraken/%s.gtdb.kraken.contig.output.gz" % ASS, nt="Analysis/taxonomy/kraken/%s.nt.kraken.contig.output.gz" % ASS
    output: temp("Analysis/taxonomy/kraken/%s.combined.kraken.contig.output.gz" % ASS)
    log: "logs/analysis/taxonomy/combine_kraken_contigs.log"
    message: "Combining GTDB and NT kraken results for contigs."

use rule make_report_base as make_report_contigs with:
    input: "Analysis/taxonomy/kraken/%s.combined.kraken.contig.output.gz" % ASS
    output: "Analysis/taxonomy/kraken/%s.kraken.contig.unfiltered.report" % ASS
    log: "logs/analysis/taxonomy/make_report_contigs.log"
    message: "Generating final combined kraken report for contigs."

use rule filter_report_base as filter_kraken_contig_report with:
    input: "Analysis/taxonomy/kraken/%s.kraken.contig.unfiltered.report" % ASS
    output: report("Analysis/taxonomy/kraken/%s.kraken.contig.report" % ASS, category="Taxonomy")
    log: "logs/analysis/taxonomy/filter_kraken_contigs_report.log"
    message: "Filtering combined Kraken report for contigs."

