====================
IMP3 output: Binning
====================




.. _output_binning_stats:

---------------------------
Stats from the binning step
---------------------------

The binning step doesn't always directly write to the ``Stats`` directory. If no taxonomy step is run, the GRiD results are
summarized together with the rest of the binning results in ``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``.
Some more results are summarized by the summary step however. 
