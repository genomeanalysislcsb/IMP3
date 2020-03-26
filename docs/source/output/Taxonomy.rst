
.. _output_taxonomy:

=====================
IMP3 output: Taxonomy
=====================

The output of the **taxonomy** step is mostly written to ``taxonomy/`` in the **Analysis** directory. The exception
is the output of GTDBtk, which is run on the bins selected by DASTool.


.. _output_taxonomy_stats:

----------------------------
Stats from the taxonomy step
----------------------------

The taxonomy step doesn't always directly write to the ``Stats`` directory. Some of the results are summarized by the summary
step, however. The GTDBtk results are summarized together with the binning results in ``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``.



