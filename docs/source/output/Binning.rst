====================
IMP3 output: Binning
====================

The IMP3 :ref:`Binning step <step_binning>` will output the results of all binners the user selected, the results from 
`DASTool <https://github.com/cmks/DAS_Tool>`_  and the results of running `GRiD <https://github.com/ohlab/GRiD>`_ on the bins that DASTool 
selected in the ``Binning`` directory within the ``outputdir`` (see :ref:`configuration <configuration>`).

-------------
Binning tools
-------------

Each binninng tool has a separate output directory within the ``Binning`` directory. Currently, IMP3 implements `MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>`_,
`MaxBin2 <https://kbase.us/applist/apps/kb_maxbin/run_maxbin2/release>`_ and :ref:`Binny <wtf_is_binny>`. Each binning tools generates a final file 
linking contigs to bins that is always named ``scaffold2bin.tsv``.

``MetaBAT`` is the output directory of `MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>`_. `MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>`_ 
outputs a tab-separated file with contig IDs and bin numbers ``metabat_res``, which is linked to in ``scaffold2bin.tsv``. In addition, the contigs in each bin, except bin *0*, are given in one *Fasta* file ``metabat_res.<bin>.fa``.

``MaxBin`` is the output directory of `MaxBin2 <https://kbase.us/applist/apps/kb_maxbin/run_maxbin2/release>`_. `MaxBin2 <https://kbase.us/applist/apps/kb_maxbin/run_maxbin2/release>`_outputs a tab-separated file with contig IDs and bin names ``maxbin_contig2bin.txt``, which is linked
to in ``scaffold2bin.tsv``. `MaxBin2 <https://kbase.us/applist/apps/kb_maxbin/run_maxbin2/release>`_ also outputs a *Fasta* file for each bin ``maxbin_res.<bin>.fasta``. Contigs which can't be put in
a bin are in ``maxbin_res.noclass`` and contigs that are too short for analysis are in ``maxbin_res.tooshort``. Summaries of the essential markers MaxBin detected in each bin and the bin size are in ``maxbin_res.marker`` and ``maxbin_res.summary``.
In addition, intermediary results are in ``maxbin_res.marker_of_each_bin.tar.gz`` and a log is recorded to ``maxbin_res.log``.

``binny`` is the output directory of :ref:`Binny <wtf_is_binny>`. :ref:`Binny <wtf_is_binny>` outputs its final result ``contigs2clusters.10.4.tsv`` and ``contigs2clusters.10.4.RDS``. The bins with at least some
completeness (names starting with P ("perfect"), G ("good"), O ("okay"), L "low completeness") are extracted and copied to ``scaffold2bin.tsv``.
:ref:`Binny <wtf_is_binny>` is based on tSNE embeddings by `VizBin <http://claczny.github.io/VizBin/>`_ of *k*-mer frequencies in contigs with masked 
(temporarily removed) rRNA genes (``<mg|mt|mgmt>.assembly.merged.cut.fa``) and the coordinates from `VizBin <http://claczny.github.io/VizBin/>`_ are stored in ``<mg|mt|mgmt>.vizbin.with-contig-names.points``. 
For comprehensibility and transparency, intermediate results are kept in ``clusterFiles.tar.gz``, ``clusterFirstScan.<pk>.<nn>.tsv``, ``bimodalClusterCutoffs.<pk>.<nn>.tsv``,
``reachabilityDistanceEstimates.<pk>.<nn>.tsv``, ``clusteringWS.<pk>.<nn>.Rdata``, and ``binny_WS.Rdata``. :ref:`Binny <wtf_is_binny>` also visualizes
its intermediary results in ``scatterPlot<1-4>.<pk>.<nn>.pdf`` and the final bins visualized in a tSNE embedding in ``finalClusterMap.10.4.png``. 

-------
DASTool
-------

The summarized output of `DASTool <https://github.com/cmks/DAS_Tool>`_ is in ``selected_DASTool_summary.txt`` and the ``selected_DASTool_bins`` directory. 
The contigs of each bin are in ``selected_DASTool_bins/<bin>.contigs.fa``. `DASTool <https://github.com/cmks/DAS_Tool>`_ uses the presence of 
single-copy marker genes to assess bins. As IMP3 has already generated gene predictions in the :ref:`Analysis step <step_analysis>`, they are used as input 
for `DASTool <https://github.com/cmks/DAS_Tool>`_ (with a changed header; ``prokka.renamed.faa``) and `DASTool <https://github.com/cmks/DAS_Tool>`_ 
keeps the annotations per gene in ``prokka.renamed.faa.<archaea|bacteria>.scg``. The lengths of all contigs are in ``selected.seqlength``.
The results of `DASTool <https://github.com/cmks/DAS_Tool>`_ assessment are in ``selected_<binner>.eval``. `DASTool <https://github.com/cmks/DAS_Tool>`_ also gives visual 
output in ``selected_DASTool_hqBins.pdf`` and ``selected_DASTool_scores.pdf`` and a log file in ``selected_DASTool.log``.

----
GRiD
----

The results of `GRiD <https://github.com/ohlab/GRiD>`_ is given for every `DASTool <https://github.com/cmks/DAS_Tool>`_ 
selected bin in ``selected_DASTool_bins/<bin>/grid``. IMP3 archives and compresses these directories after 
adding the information into ``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``. All the archived data is in ``per_bin_results.tar.gz``.

------
GTDBtk
------

If the :ref:`Taxonomy step <step_taxonomy>` is selected in addition to the :ref:`Binning step <step_binning>`, the selected bins from 
`DASTool <https://github.com/cmks/DAS_Tool>`_ are analysed with `GTDBtk <https://github.com/Ecogenomics/GTDBTk>`_. The results are 
in ``selected_DASTool_bins/<bin>/GTDB``. IMP3 archives this data in ``per_bin_results.tar.gz``, after the results are added to
``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``.

.. _output_binning_stats:

---------------------------
Stats from the Binning step
---------------------------

The :ref:`Binning step <step_binning>` writes to the ``Stats`` directory. If no :ref:`Taxonomy step <step_taxonomy>` is run, the 
`GRiD <https://github.com/ohlab/GRiD>`_ results are summarized together with the rest of the binning results in ``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``.
If taxonomy is performed, the `GRiD <https://github.com/ohlab/GRiD>`_ results are also combined with the `GTDBtk <https://github.com/Ecogenomics/GTDBTk>`_ 
results in ``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``. Some more results are summarized by the :ref:`Summary step <step_summmary>` if defined. 
