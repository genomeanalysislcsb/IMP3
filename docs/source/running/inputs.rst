.. _input_options:

==================
IMP3 input options
==================

IMP3 is designed to perform integrated analyses of **metaG** and **metaT** data. Inputs can be either raw sequencing reads or already assembled contigs 
and mapped reads.

A typical :ref:`workflow <steps_overview>` starts with a pair of files of paired-end **metaG** reads and a pair
of **metaT** reads (both in `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`_ format, either gzipped or not). Alternatively, IMP3 can take only **metaG**
or only **metaT** reads. 

If the data is already assembled and the contigs should be annotated and binned into metagenomics-assembled genomes (**MAGs**), 
IMP3 also takes assemblies (in FASTA format) in addition to alignments (`BAM <https://genome.sph.umich.edu/wiki/BAM>`_ files) or trimmed reads (`FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`_ format) 
as input.

All inputs are defined in the :ref:`config file <configuration>`.

-----
Reads
-----

The ``Metagenomics`` and ``Metatranscriptomics`` input fields expect two or three `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`_ or
gzipped `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`_ files, 
separated by a space. The first two files should be the forward and reverse reads, if IMP3 analysis should start from pre-processing reads. The reads 
need to be in the same order in both files. 

For processing the original raw `Fastq <https://en.wikipedia.org/wiki/FASTQ_format>`_ files, the following setting should be used:

.. code-block:: yaml

  raws: 
    Metagenomics: "/path/to/metagenomics.read.r1.fastq.gz /path/to/metagenomics.read.r2.fastq.gz"
    Metatranscriptomics: "/path/to/metatranscriptomics.read.r1.fastq.gz /path/to/metatranscriptomics.read.r2.fastq.gz"
    LongReads: ""
    LongReadTech: ""
    Contigs: ""
    Alignment_metagenomics: ""
    Alignment_metatranscriptomics: ""

If the user has **already pre-processed** reads, either for the :ref:`Assembly <step_assembly>` step or for further analysis,
the user is required to pass THREE files, two paired read files and one additional single ends file. The singleton read file is given last.

**Note**: *In this case, the user should NOT include* ``preprocessing`` *in the IMP* :ref:`steps <steps_overview>`.

.. code-block:: yaml

  raws: 
    Metagenomics: "/path/to/processed_metagenomics.read.r1.fastq.gz /path/to/processed_metagenomics.read.r2.fastq.gz /path/to/processed_metagenomics.read.se.fastq.gz"
    Metatranscriptomics: "/path/to/processed_metatranscriptomics.read.r1.fastq.gz /path/to/processed_metatranscriptomics.read.r2.fastq.gz /path/to/processed_metatranscriptomics.read.se.fastq.gz"
    LongReads: ""
    LongReadTech: ""
    Contigs: ""
    Alignment_metagenomics: ""
    Alignment_metatranscriptomics: ""
    
The user can also give only **metaG** or only **metaT** reads to IMP3, either raw or pre-processed.

If the user wishes to perform a **metaG** assembly including long reads, a single file of **already pre-processed** long-reads data (in `Fastq <https://en.wikipedia.org/wiki/FASTQ_format>`_ format). In this case, additionally the sequencing method 
(possible values are ``nanopore`` and ``pacbio``) can be added.

**Note**: *The long reads are only used by IMP3, if* ``metaspades`` *is chosen as assembler.*

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

If the user has already an assembly and would like to use IMP3 to annotate genes, perform binning and/or determine contig-level
taxonomy, the contigs can be used as input in FASTA format (NOTE: the FASTA header should only contain the contig name, NO other information and NO spaces: >CONTIGNAME1). 
In addition, the user can give reads either in `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`_ files or already aligned
as `BAM <https://genome.sph.umich.edu/wiki/BAM>`_ files. For `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`_ files, the same limitations apply as discussed above. 
The `BAM <https://genome.sph.umich.edu/wiki/BAM>`_ files should be sorted by contig name and coordinate, and again, the contig names are NOT allowed to contain spaces !!!).

.. code-block:: yaml

  raws: 
    Metagenomics: ""
    Metatranscriptomics: ""
    LongReads: ""
    LongReadTech: ""
    Contigs: "/path/to/assembly.fasta"
    Alignment_metagenomics: "/path/to/metagenomics/read.sorted.alignment.bam"
    Alignment_metatranscriptomics: "/path/to/metatranscriptomics/read.sorted.alignment.bam"


