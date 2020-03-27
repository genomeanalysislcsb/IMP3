====================
IMP3 output: Binning
====================

IMP3 will output the results of all the binners you selected, what `DASTool <https://github.com/cmks/DAS_Tool>`_ made of the bins and
the results of running `GRiD <https://github.com/ohlab/GRiD>`_ on the bins that DASTool selected in the **Binning** directory.

-------
Binners
-------

Each binner has one output directory. Currently, IMP3 implements `MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>`_,
`MaxBin2 <https://kbase.us/applist/apps/kb_maxbin/run_maxbin2/release>`_ and :ref:`binny <wtf_is_binny>`. The final file
linking contigs to bins is always named ``scaffold2bin.tsv``.

**MetaBAT** simply outputs a tab-separated file with contig IDs and bin numbers ``metabat_res``, which is linked to in ``scaffold2bin.tsv``.
In addition, the contigs in each bin, except 0, are given in one FASTA file ``metabat_res.<bin>.fa``.

**MaxBin** similarly outputs a tab-separated file with contig IDs and bin names ``maxbin_contig2bin.txt``, which is linked
to in ``scaffold2bin.tsv``. MaxBin also outputs a FASTA file for each bin ``maxbin_res.<bin>.fasta``. Contigs which can't be put in
a bin are in ``maxbin_res.noclass`` and contigs that are too short for analysis are in ``maxbin_res.tooshort``.
Summaries of the essential markers MaxBin detected in each bin and the bin size are in ``maxbin_res.marker`` and ``maxbin_res.summary``.
In addition, intermediary results are in ``maxbin_res.marker_of_each_bin.tar.gz`` and a log is recorded to ``maxbin_res.log``.

**binny** outputs its final result ``contigs2clusters.10.4.tsv`` and ``contigs2clusters.10.4.RDS``. The bins with at least some
completeness (names starting with P ("perfect"), G ("good"), O ("okay"), L "low completeness") are extracted and copied to ``scaffold2bin.tsv``.
Binny is based on tSNE embeddings by `VizBin <http://claczny.github.io/VizBin/>`_
of *k*-mer frequencies in contigs with masked (temporarily removed) rRNA genes (``<mg|mt|mgmt>.assembly.merged.cut.fa``)
and the coordinates from VizBin are stored in ``<mg|mt|mgmt>.vizbin.with-contig-names.points``. 
For comprehensibility and transparency, intermediate results
are kept in ``clusterFiles.tar.gz``, ``clusterFirstScan.<pk>.<nn>.tsv``, ``bimodalClusterCutoffs.<pk>.<nn>.tsv``,
``reachabilityDistanceEstimates.<pk>.<nn>.tsv``, ``clusteringWS.<pk>.<nn>.Rdata``, and ``binny_WS.Rdata``. Binny also visualizes
its intermediary results in ``scatterPlot<1-4>.<pk>.<nn>.pdf`` and the final bins visualized in a tSNE embedding in
``finalClusterMap.10.4.png``. 

-------
DASTool
-------

The summarized output of `DASTool <https://github.com/cmks/DAS_Tool>`_ is in ``selected_DASTool_summary.txt``.
The contigs of each bin are in ``selected_DASTool_bins/<bin>.contigs.fa``. DASTool uses the presence of single-copy
marker genes to assess bins. As IMP3 has already produced gene predictions, we give these to DASTool (with a changed header;
``prokka.renamed.faa``) and DASTool keeps the annotations per gene in ``prokka.renamed.faa.<archaea|bacteria>.scg``.
The lengths of all contigs are in ``selected.seqlength``.
The results of DASTools assessment are in ``selected_<binner>.eval``. DASTool also gives visual output in ``selected_DASTool_hqBins.pdf``
and ``selected_DASTool_scores.pdf`` and a log file in ``selected_DASTool.log``.

----
GRiD
----

The result of **GRiD** is produced in ``selected_DASTool_bins/<bin>/grid``. However, IMP3 archives and compresses this directory
after writing the information into ``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``. All the archived data is in ``per_bin_results.tar.gz``.

------
GTDBtk
------

If the **taxonomy** step is selected in addition to binning, the bins from DASTool are analysed
with `GTDBtk <https://github.com/Ecogenomics/GTDBTk>`_. The results are produced in ``selected_DASTool_bins/<bin>/GTDB``, but
just like with GRiD, the data is archived in ``per_bin_results.tar.gz``, after the results are recorded in
``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``


.. _output_binning_stats:

---------------------------
Stats from the binning step
---------------------------

The binning step writes to the ``Stats`` directory. If no taxonomy step is run, the GRiD results are
summarized together with the rest of the binning results in ``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``.
If taxonomy is performed, the GRiD results are also combined with the GTDBtk results in
``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``.
Some more results are summarized by the summary step. 
