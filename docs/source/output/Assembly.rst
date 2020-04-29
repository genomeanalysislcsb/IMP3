.. _output_assembly:

=====================
IMP3 output: Assembly
=====================

During the :ref:`Assembly step <step_assembly>` (or if the user provided an existing assembly, see :ref:`input <input_options>`) all outputs will be written to 
the ``Assembly`` directory within the defined ``outputdir`` directory (see :ref:`configuration <configuration>`). 

- The final output is a FASTA file of the assembled contigs ``<mg|mt|mgmt>.assembly.merged.fa`` [\*] . The FASTA headers contain the **sample name** as given in the :ref:`config file <configuration>`, separated by an underscore, *contig*, another underscore, and a number, e.g. *test_contig_1*. 
- The **index** files of the final contig file will be generated , namely for `BWA <http://bio-bwa.sourceforge.net/>`_ (suffixes ``amb``, ``ann``, ``bwt``, ``pac``, and ``sa``), `Samtools <http://www.htslib.org/doc/faidx.html>`_ (``fai``) and bioperl (``index``). 
- A ``bed3`` file is also stored for later access.
**Note**: *some of these files are produced by the* :ref:`Analysis step <step_analysis>` *, so they will not be present after running only the* :ref:`Assembly step <step_assembly>`.
 
During the :ref:`Assembly step <step_assembly>`, IMP3 maps back the **processed** reads to the final contigs and stores the alignment 
``<mg|mt>.reads.sorted.bam`` and index ``<mg|mt>.reads.sorted.bam.bai``. The `BAM <https://genome.sph.umich.edu/wiki/BAM>`_ files are sorted
by contig name and position.
 
The :ref:`Assembly step <step_assembly>` actually consists of a large number of sub-steps (:ref:`iterative assembly <step_assembly>`) generating a huge 
amount of intermediate result files and directories that will be archived and compressed into ``intermediary.tar.gz``.

.. _output_assembly_stats:

------------------------
Stats from Assembly step
------------------------

The :ref:`Assembly step <step_assembly>` will also collect some summary statistics and save them to the ``Stats`` directory: 

- The GC-content of the final contigs is recorded in the tab-separated file ``<mg|mt|mgmt>/<mg|mt|mgmt>.assembly.gc_content.txt``, which contains a header and holds two columns for the contig names and GC in percent (i.e 0-100), respectively. 
- The length of the contigs is provided in a tab-separated file ``<mg|mt|mgmt>/<mg|mt|mgmt>.assembly.length.txt`` with two simple columns, contig names and lengths.
- The stats on the read mapping are kept in the ``Stats`` subdirectories with the **metaG** and/or **metaT** (``mg/`` or ``mt/``):
- ``<mg|mt>/<mg|mt|mgmt>.assembly.contig_flagstat.txt`` contains the numeric part of the `samtools flagstat <http://www.htslib.org/doc/samtools-flagstat.html>`_ output. 
- the average depth of coverage for each set of reads for all contigs that have at least one read mapping to them is given in ``<mg|mt>/<mg|mt|mgmt>.assembly.contig_depth.txt``. This file is headerless and tab-separated, with the contig names in the first column and the average depth of coverage in the second. 
- The :ref:`Binning step <step_binning>` adds based on this the file ``<mg|mt>/<mg|mt|mgmt>.assembly.contig_depth.0.txt``, which also contains lines for contigs with zero coverage. 
- The file ``<mg|mt>/<mg|mt|mgmt>.assembly.contig_coverage.txt`` contains the processed output of `bedtools genomeCoverageBed <https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html>`_. It contains seven tab-separated columns: 
    - 1: contig names,
    - 2: 0 (= the first position in bed coordinates),
    - 3: contig length (= the last position in bed coordinates),
    - 4: the number of regions overlapping with the aligned positions in the bam file (roughly = number of mapping reads),
    - 5: the number of covered positions, 
    - 6: the number of positions (= value in column 3),
    - 7: the coverage breadth, i.e. proportion of covered positions (scaling 0-1)).


.. [*] The name of the assembly depends on the workflow: if the workflow has only **metaG** reads as input, the user has chosen the  **non-hybrid** assembly 
   workflow, or has provided an exsisting assembly, the assembly will be referred to by ``mg``. If only **metaT** reads were given, the assembly will be 
   referred to by ``mt``. If the **hybrid** **metaG** and **metaT** workflow was defined and both types of reads were supplied, the assembly will be
   represented by ``mgmt``.

