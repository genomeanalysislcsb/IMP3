
.. _output_preprocessing:

==========================
IMP3 output: Preprocessing
==========================

The :ref:`preprocessing step <step_preprocessing>` creates standard .fastq files in the ``Preprocessing`` directory and also
in ``Stats``. The exact identity of the files depends on your settings. 

In all cases, **the final set of reads** 
(that will go into the :ref:`assembly step <step_assembly>`) will be pointed to by symbolic
links ``<mg|mt>.<r1|r2|se>.preprocessed.fq`` [#]_ , [#]_. Your input reads will be also be in the directory
(called ``<mg|mt>.<r1|r2>.fq.gz``). All IMP3 runs with preprocessing will also have the
`trimmomatic output <http://www.usadellab.org/cms/?page=trimmomatic>`_ ``<mg|mt>.<r1|r2|se>.trimmed.fq.gz``. 

Metatranscriptomic reads are always filtered to remove rRNA reads using `SortMeRNA <https://github.com/biocore/sortmerna>`_ .
The remaining reads from this step are ``mt.<r1|r2|se>.trimmed.rna.fq.gz``.
If you chose to filter your reads by mapping against a reference, you will also find ``mt.<r1|r2|se>.trimmed.rna_filtered.<reference>_filtered.fq.gz``
and/or ``mg.<r1|r2|se>.trimmed.<reference>_filtered.fq.gz``.

**Why is there a Preprocessing directory if you did not choose to do preprocessing?** If you supplied already preprocessed reads,
these will be placed in the Preprocessing directory, because that's where the next steps (assembly or analysis) will look for them.


.. _output_preprocessing_stats:

-----------------------------
Stats from preprocessing step
-----------------------------

In addition to the preprocessing, the preprocessing step will also collect some stats for you and save them to the ``Stats``
directory. To this end, it runs `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ on the input reads and preprocessed reads.
The resulting .html report and directories with more detailed tables are in the respective
subdirectories ``<mg|mt>/`` named ``<mg|mt>.<r1|r2>_fastqc.<zip|html>`` and ``<mg|mt>.<r1|r2|se>.preprocessed_fastqc.<zip|html>``.


.. [#] ``mg`` denotes metagenomics data throughout the workflow, ``mt`` represents metatranscriptomics data.
.. [#] ``r1`` contains first reads with a partner, which will be in ``r2``, ``se`` contains the single ends that have lost their partner to preprocessing.




