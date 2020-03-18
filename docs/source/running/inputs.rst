=================
IMP input options
=================

IMP is designed to perform integrated analyses of metagenomics and metatranscriptomics data.
A typical workflow starts with a pair of files of paired-end metagenomics reads and a pair
of files with paired-end metatranscriptomics reads (both in fastq format). Alternatively, you can hand IMP only metagenomics
or only metatranscriptomics reads. If you want to forego assembly and are looking to annotate
your assembly and to bin contigs into metagenomics-assembled genomes, you can give assemblies (as fasta)
in addition to alignments (.bam files) or trimmed reads (fastq format) to IMP.

All inputs are defined in the :ref:`config file <configuration>`. The Metagenomics and Metatranscriptomics fields expect
two or three fastq or gzipped fastq files, separated by a space. These files should be the forward and reverse reads,
if you start the analysis from pre-processing reads. If you supply already pre-processed reads, either to the
assembly step or for further analysis, you can pass one additional single ends file.

.. code-block:: yaml

raws: 
  Metagenomics: "/path/to/metagenomics.read.r1.fastq.gz /path/to/metagenomics.read.r2.fastq.gz"
  Metatranscriptomics: "/another/path/to/metatranscriptomics.read.r1.fastq.gz /another/path/to/metatranscriptomics.read.r2.fastq.gz"
  LongReads: ""
  LongReadTech: ""
  Contigs: ""
  Alignment_metagenomics: ""
  Alignment_metatranscriptomics: ""


If you wish to perform a metagenomics assembly, you can also add a single file of already pre-processed
long-reads data (as fastq files). In this case, please also name the sequencing method (possible values are "nanopore" and "pacbio").

.. code-block:: yaml

raws: 
  Metagenomics: "/path/to/metagenomics.read.r1.fastq.gz /path/to/metagenomics.read.r2.fastq.gz"
  Metatranscriptomics: ""
  LongReads: "/path/to/longread/reads.fastq"
  LongReadTech: ""
  Contigs: ""
  Alignment_metagenomics: ""
  Alignment_metatranscriptomics: ""

 
