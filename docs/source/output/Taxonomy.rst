.. _output_taxonomy:

=====================
IMP3 output: Taxonomy
=====================

The output of the :ref:`taxonomy step <step_taxonomy>` is mostly written to ``taxonomy`` in the **Analysis** directory. The exception
is the output of `GTDBtk <https://github.com/Ecogenomics/GTDBTk>`_, which is run on the bins selected by `DASTool <https://github.com/cmks/DAS_Tool>`_ 
and has its results summarized with the bins.


.. _output_taxonomy_stats:

----------------------------
Stats from the taxonomy step
----------------------------

The :ref:`Taxonomy step <step_taxonomy>` does not always directly write to the ``Stats`` directory. Some of the results are summarized by the :ref:`Summary <step_summary>`. 
The `GTDBtk <https://github.com/Ecogenomics/GTDBTk>`_ results are summarized together with the results from the :ref:`binning step <step_binning>` results in ``Stats/<mg|mt|mgmt>/<mg|mt|mgmt>.bins.tsv``.



