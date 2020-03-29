
.. _output_assembly:

=====================
IMP3 output: Assembly
=====================

If you run the assembly step (or provide an assembly as input), a FASTA file of the contigs will be kept in **Assembly** named
``<mg|mt|mgmt>.assembly.merged.fa`` [*]_ . The headers contain the sample name, separated by an underscore from "contig", another
underscore, and a number. In addition, a few indices of the contigs will be stored here, namely for
`BWA <http://bio-bwa.sourceforge.net/>`_ (suffixes ``amb``, ``ann``, ``bwt``, ``pac``, and ``sa``),
`Samtools <http://www.htslib.org/doc/faidx.html>`_ (``fai``) and bioperl (``index``). A ``bed3`` file is also stored for later access.
*Note: some of these files are produced by the analysis step, so you may not see all of them immediately or ever.*
 
Next to the assembly, the assembly step will map the reads back to the contig and keep the alignment ``<mg|mt>.reads.sorted.bam``
and index ``<mg|mt>.reads.sorted.bam.bai``. The reads are sorted by coordinate.
 
The assembly step actually consists of a large number of steps (:ref:`iterative assembly <step_assembly>`). You will see
a number of files being written into the directory if you watch IMP3 work. These intermediary files are kept,
but archived and compressed in ``intermediary.tar.gz``.


.. _output_assembly_stats:

------------------------
Stats from assembly step
------------------------

In addition to running the assembly and mapping, the assembly step will also collect some stats for you and save them to the **Stats**
directory. In particular, GC-content is recorded in the tab-separated file ``<mg|mt|mgmt>/<mg|mt|mgmt>.assembly.gc_content.txt``,
which contains a header and holds two columns for the contig names and GC in percent (i.e 0-100), respectively. And the
length of the contigs is provided in a tab-separated file ``<mg|mt|mgmt>/<mg|mt|mgmt>.assembly.length.txt`` with two simple columns,
contig names and lengths.

The stats of the read mapping are kept in the subdirectory with the metagenomic and/or metatranscriptomic data (``mg/`` or ``mt/``).
``<mg|mt>/<mg|mt|mgmt>.assembly.contig_flagstat.txt`` contains the numeric part of
the `samtools flagstat <http://www.htslib.org/doc/samtools-flagstat.html>`_ output. The average depth of coverage with each
set of reads is recorded for all contigs that have at least one read mapping to them in ``<mg|mt>/<mg|mt|mgmt>.assembly.contig_depth.txt``. This
file is headerless and tab-separated, with the contig names in the first column and the average depth of coverage in the second.
(The binning step adds based on this the file ``<mg|mt>/<mg|mt|mgmt>.assembly.contig_depth.0.txt`` which also contains lines for contigs with
0 coverage). Finally, the file ``<mg|mt>/<mg|mt|mgmt>.assembly.contig_coverage.txt`` is the output of
`bedtools' genomeCoverageBed <https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html>`_ (plus a little AWK). It contains seven tab-separated
columns (1: contig names, 2: 0 (= the first position in bed coordinates), 3: contig length (= the last position in bed coordinates),
4: the number of regions overlapping with the aligned positions in the bam file (roughly = number of mapping reads),
5: the number of covered positions, 6: the number of positions (= value in column 3),
7: the coverage breadth, i.e. proportion of covered positions (scaling 0-1)).


.. [*] The name of the assembly depends on the workflow: if you have only metagenomic reads, have chosen the non-hybrid assembly workflow,
   or have provided your own assembly, the assembly will be referred to by ``mg``. If you used only metatranscriptomic reads, it will be 
   referred to by ``mt``. If you chose the hybrid metagenomics and metatranscriptomics workflow and supplied both kinds of reads, it will be
   represented by ``mgmt``.

