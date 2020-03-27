====================
IMP3 output: Binning
====================

DASTool result:	selected_ DASTool_summary.txt
contigs per bin:	selected_ DASTool_bins/<bin>.contigs.fa
protein predictions:	 prokka.renamed.faa
single-copy ann.:	prokka.renamed.faa.<taxon>.scg
binner inputs:	selected_<binner>.eval
others:	selected_ DASTool_<data>.<suffix>
Binner outputs:	<binner>/
Binner final results:	<binner>/scaffold2bin.tsv
More results per bin:	per_bin_results.tar.gz



.. _output_binning_stats:

---------------------------
Stats from the binning step
---------------------------

The binning step doesn't always directly write to the ``Stats`` directory. If no taxonomy step is run, the GRiD results are
summarized together with the rest of the binning results in ``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``.
Some more results are summarized by the summary step however. 
