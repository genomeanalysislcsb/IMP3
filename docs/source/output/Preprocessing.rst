.. _output_preprocessing:

==========================
IMP3 output: Preprocessing
==========================

The :ref:`Preprocessing step <step_preprocessing>` creates output files in the ``Preprocessing`` and ``Stats`` directories within the defined ``outputdir``. 
The exact naming of the files depends on the :ref:`configuration <configuration>` settings. 

For all runs with reads as :ref:`input <input_options>`, the original input reads will be copied to the ``Preprocessing`` 
directory and renamed into ``<mg|mt>.<r1|r2>.fq.gz``. **The final set of processed reads** (that will be the input for the :ref:`Assembly step <step_assembly>`) 
will be pointed to by symbolic links ``<mg|mt>.<r1|r2|se>.preprocessed.fq`` [*]_ , [*]_. All IMP3 runs with the step :ref:`Preprocessing <step_preprocessing>` 
will also have the trimmed reads processed with `trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_ (``<mg|mt>.<r1|r2|se>.trimmed.fq.gz``). 

**MetaT** reads are always filtered to remove rRNA reads using `SortMeRNA <https://github.com/biocore/sortmerna>`_ . The remaining reads from this filtering step 
are named ``mt.<r1|r2|se>.trimmed.rna.fq.gz``. If reads were mapped against one or more reference genomes, the filtered reads are named 
``mt.<r1|r2|se>.trimmed.rna_filtered.<reference>_filtered.fq.gz`` and/or ``mg.<r1|r2|se>.trimmed.<reference>_filtered.fq.gz``.

**Note**: *Why is there a Preprocessing directory if you did not choose to do preprocessing?* If already preprocessed reads were given as :ref:`input <input_options>`,
these reads will be placed in the ``Preprocessing`` directory following the exact same naming conventions as described above, to ensure that the downstream 
steps can use the same consistent input files and directories.


.. _output_preprocessing_stats:

-----------------------------
Stats from Preprocessing step
-----------------------------

The :ref:`Preprocessing step <step_preprocessing>` will also collect some statistics on the original and processed reads and save them to the **Stats**
directory. To this end, it runs `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ on the original and preprocessed reads. The resulting HTML reports 
and directories with more detailed tables are in the respective subdirectories ``<mg|mt>/`` named ``<mg|mt>.<r1|r2>_fastqc.<zip|html>`` and ``<mg|mt>.<r1|r2|se>.preprocessed_fastqc.<zip|html>``.
In addition, the number of reads in each step is counted in ``<mg|mt>/<mg|mt>.read_counts.txt``, which holds two tab-separated columns with the
file names and number of reads, respectively.


.. [*] ``mg`` denotes **metaG** data throughout the workflow, ``mt`` represents **metaT** data.
.. [*] ``r1`` contains first reads with a partner, which will be in ``r2``, ``se`` contains the single ends that have lost their partner during the :ref:`Preprocessing step <step_preprocessing>`.




