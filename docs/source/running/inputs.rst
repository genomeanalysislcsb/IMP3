.. _input_options:

==================
IMP3 input options
==================

IMP3 is designed to perform integrated analyses of metagenomics and metatranscriptomics data.
A typical workflow starts with a pair of files of paired-end metagenomics reads and a pair
of files with paired-end metatranscriptomics reads (both in fastq format). Alternatively, you can hand IMP3 only metagenomics
or only metatranscriptomics reads. If you want to forego assembly and are looking to annotate
your assembly and to bin contigs into metagenomics-assembled genomes, you can give assemblies (as fasta)
in addition to alignments (.bam files) or trimmed reads (fastq format) to IMP3.

All inputs are defined in the :ref:`config file <configuration>`.

-----
Reads
-----

The Metagenomics and Metatranscriptomics fields expect
two or three fastq or gzipped fastq files, separated by a space. These files should be the forward and reverse reads,
if you start the analysis from pre-processing reads. The read partners need to be in the same order in both files, singletons are not
permitted.

.. code-block:: yaml

  raws: 
    Metagenomics: "/path/to/metagenomics.read.r1.fastq.gz /path/to/metagenomics.read.r2.fastq.gz"
    Metatranscriptomics: "/another/path/to/metatranscriptomics.read.r1.fastq.gz /another/path/to/metatranscriptomics.read.r2.fastq.gz"
    LongReads: ""
    LongReadTech: ""
    Contigs: ""
    Alignment_metagenomics: ""
    Alignment_metatranscriptomics: ""

If you supply already pre-processed reads, either to the assembly step or for further analysis,
you can pass one additional single ends file. *In this case, you should not include "preprocessing" in the IMP steps.*

.. code-block:: yaml

  raws: 
    Metagenomics: "/path/to/processed_metagenomics.read.r1.fastq.gz /path/to/processed_metagenomics.read.r2.fastq.gz /path/to/processed_metagenomics.read.se.fastq.gz"
    Metatranscriptomics: "/another/path/to/processed_metatranscriptomics.read.r1.fastq.gz /another/path/to/processed_metatranscriptomics.read.r2.fastq.gz /another/path/to/processed_metatranscriptomics.read.se.fastq.gz"
    LongReads: ""
    LongReadTech: ""
    Contigs: ""
    Alignment_metagenomics: ""
    Alignment_metatranscriptomics: ""
    
You can also give only metagenomics or only metatranscriptomics reads to IMP3.

If you wish to perform a metagenomics assembly including long reads, you can also add a single file of **already pre-processed**
long-reads data (as fastq files). In this case, please also name the sequencing method (possible values are "nanopore" and "pacbio").
*The long reads are only used by IMP3, if you choose "metaspades" as assembler.*

.. code-block:: yaml

  raws: 
    Metagenomics: "/path/to/metagenomics.read.r1.fastq.gz /path/to/metagenomics.read.r2.fastq.gz"
    Metatranscriptomics: ""
    LongReads: "/path/to/longread/reads.fastq"
    LongReadTech: ""
    Contigs: ""
    Alignment_metagenomics: ""
    Alignment_metatranscriptomics: ""

----------------------
Contigs and alignments
----------------------

If you already have an assembly and would like to use IMP3 to annotate genes, perform binning and/or determine contig-level
taxonomy, you can input the contigs in .fasta format. In addition, you can give reads either as .fastq files or already aligned
as .bam files. For .fastq-files, the same limitations apply as discussed above. The alignment should be sorted by name.

.. code-block:: yaml

  raws: 
    Metagenomics: ""
    Metatranscriptomics: ""
    LongReads: ""
    LongReadTech: ""
    Contigs: "/path/to/assembly.fasta"
    Alignment_metagenomics: "/path/to/metagenomics/read.sorted.alignment.bam"
    Alignment_metatranscriptomics: "/path/to/metatranscriptomics/read.sorted.alignment.bam"


