.. _configuration:

=====================
Configuring IMP3 runs
=====================

IMP3 has a modular design and is highly configurable. Users can define the analysis steps IMP3 should perform, the data it should run on,
the databases and tools that should be used and the settings of each step in a single **config file** in `yaml <https://yaml.org/>`_ format.
This **config file** is provided to `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ 
as described :ref:`here <run_IMP>`, using the ``--configfile`` argument.

This description guides through all sections and fields of the **config file**. 

Some general notes on the **config file**:

- The **order of the items** in the config file is **arbitrary**.
- The config file does **not** need to contain all fields:
- only two fields are **required**: :ref:`raws <input_options>` and an **output directory**. 
- all other fields are optional and if emtpy the IMP3 default settings will be used. 
- IMP3 will write the user-specified configuration including pre-defined fields to the output folder as ``sample.config.yaml``.

-------------------------------
Inputs, outputs and directories
-------------------------------

The user has to tell IMP3 where the **raw or input files** are located. (The default is no input, and nothing will be done). 
Different :ref:`steps <steps_overview>` in IMP3 require different input files. All input files are given as absolute paths.
Input files can be:

- ``Metagenomics`` and/or ``Metatranscriptomics`` paired-end or single-end raw `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`_ files. Paired end read files are separated by ``<space>``, or
- ``Alignment_metagenomics`` and/or ``Alignment_metatranscriptomics``: MetaG and/or MetaT `BAM <https://genome.sph.umich.edu/wiki/BAM>`_ files (requires ``Contigs`` file, required if the :ref:`Preprocessing <step_preprocessing>` and :ref:`Assembly <step_assembly>` steps are skipped).
- ``Contigs`` file in FASTA format (used for the mapping of the metaG und metaG raw files, required if the :ref:`Preprocessing <step_preprocessing>` and :ref:`Assembly <step_assembly>` steps are skipped).
- ``LongReads`` file in `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`_ (requires a ``LongReadTech`` description). ``LongReads`` are currently only supported with additional short **metaG** reads and with Spades assembler.
- ``LongReadTech`` description of used long read sequencing technology, currently ``nanopore`` or ``pacbio``.
- ``Gff`` file based on the ``Contigs`` file (only required if the :ref:`Preprocessing <step_preprocessing>`, :ref:`Assembly <step_assembly>` and :ref:`Analysis <step_analysis>` steps are skipped and only the step :ref:`Binning <step_binning>` should performed).

The ``outputdir`` (with full path) is a required field and defines the directory the IMP3 workflow output should be written to.

Users can give a **sample name** (``sample``, default value is ``test``). Special characters should be avoided, punctuation will be converted to underscore.
The **sample name** will be used throughout the IMP3 workflow, e.g. prepended to contig names. If the user-defined **sample name** is empty ``sample: ""``, IMP3 
will take the last two parts of the output path, concatenate them with an underscore and use the results as sample name.

A directory for **temporary files** can be given (``tmp_dir``, default name is ``tmp``). If an absolute path is not given, the ``tmp`` directory will be placed into the output 
directory.

If IMP3 performs any of the steps that need :ref:`databases <ext_DBs>` (i.e. :ref:`Preprocessing <step_preprocessing>` (filtering of reads against a host genome, filtering rRNA reads),
:ref:`Analysis <step_analysis>` (gene calling, functional annotation), or :ref:`Taxonomy <step_taxonomy>` (kraken2), the user has to give the full path to a single directory with all the
databases.

.. code-block:: yaml

  raws: 
    Metagenomics: "/path/to/pair1 /path/to/pair2"
    Metatranscriptomics: "/path/to/pair1 /path/to/pair2"
    LongReads: "/path/to/longreads"
    LongReadTech: "nanopore|pacbio"
    Contigs: "/path/to/contigfile"
    Alignment_metagenomics: "/path/to/mg-bamfile"
    Alignment_metatranscriptomics: "/path/to/mt-bamfile"
    Gff: "/path/to/gfffile"
    
  outputdir: "/path/to/output"
  sample: <sampleid> (default: test)
  tmp_dir: /path/to/tmp (default: tmp)
  db_path: "/path/to/IMP_DBs"


----------
IMP3 steps
----------

IMP3 runs either the whole workflow or just some :ref:`steps <steps_overview>`. The field ``steps`` lists the user-defined IMP3 analysis steps, separated by
a space. The default is to run all steps.

The first step in IMP3 is :ref:`Preprocessing <step_preprocessing>`, i.e. quality of reads, trimming and the removal of reads which map against a reference genome, 
most commonly a host genome). The default setting (``preprocessing_filter: true``), given that **metaG** or **metaT** reads are provided, will quality-trim reads and 
map them mapped against the human genome given as file ``hg38.fa`` in the database directory under ``db_path``. The user can set ``preprocessing_filtering`` to 
``false`` to skip this step, and the trimmed, unfiltered reads will directly go into the :ref:`Assembly <step_assembly>` step.

If ``steps`` contains ``summary``, IMP3 will summarize the results. The user can define which steps should be taken for summary in the field ``summary_steps``. 
The default is to extract summary statistics (``stats``) on all performed steps and to visualize them (``vis``). If ``summary_steps`` is empty, IMP3 will 
not make summaries. Currently, the user can either choose ``"stats"`` or ``"stats vis"``. 

.. code-block:: yaml

  steps: "preprocessing assembly analysis binning taxonomy summary"
  preprocessing_filtering: true
  summary_steps: "stats vis"

.. _hardware:

------------------
Available hardware
------------------

Some IMP3 steps need high memory cpus depending on the size of the dataset and the databases used. IMP3 is able to run tasks on specific 
type of nodes. High memory demanding jobs are automatically scheduled to ``big_mem`` nodes. 
Depending on how much **RAM** on the local computer or compute cluster are available for the IMP3 run, the user can define the settings for  **normal** 
and **bigmem** cpus:

- ``normal_mem_per_core_gb`` max. available RAM (in Gb) per core of **normal** compute cores, default 4 Gb
- ``big_mem_total_gb`` max. available total RAM (in Gb) of **bigmem** compute cores, default 1600 Gb
- ``big_mem_cores`` available number of **bigmem** compute cores, default 8
- ``big_mem_per_core`` max. available RAM (in Gb) per core of **bigmem** compute cores, default 200 Gb

The number of cores that are actually used at any point in time is determined by the ``--cores`` argument of the `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ :ref:`call <run_IMP>`, 
rather than set here. 

**Note**: *Defaults may be blatantly inappropriate for your system. Most steps will run with parameters that don't match
your system, but some may cause errors (e.g. mapping, kraken, assembly).*

Currently you need to have both kinds of cores, if you run IMP3 in :ref:`cluster mode <run_cluster>`, or you have to change the :ref:`cluster config file <cluster_config>`. 

.. code-block:: yaml

  mem:
    normal_mem_per_core_gb: 4
    big_mem_total_gb: 1600
    big_mem_cores: 8
    big_mem_per_core_gb: 200


------------------------
Pre-processing: trimming 
------------------------

For the IMP3 :ref:`Preprocessing <step_preprocessing>` step, IMP3 uses `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_ for read trimming.
The user can set all ``trimmomatic`` options (defaults given below).

If the reads originate from an `Illumina Nextseq <https://www.illumina.com/systems/sequencing-platforms/nextseq.html?langsel=/us/>`_ machine, we advise 
to set ``nextseq`` to true. This will remove trailing G's that the `Nextseq <https://www.illumina.com/systems/sequencing-platforms/nextseq.html?langsel=/us/>`_ systems 
add to reads that are shorter than the run cycles.

.. code-block:: yaml

  trimmomatic:
    adapter: "TruSeq3-PE"
    leading: 20
    minlen: 40
    palindrome_clip_threshold: 30
    simple_clip_threshold: 10
    trailing: 20
    seed_mismatch: 2
    window_size: 1
    window_quality: 3
    strictness: 0.5
    target_length: 40
  nextseq: false


-------------------------
Pre-processing: filtering 
-------------------------

The IMP3 :ref:`Preprocessing <step_preprocessing>` step allows to map against a user-defined reference genome (e.g. host). 
The user can set the file name (minus suffix) of the reference genome (``filter``), which should be located or (soft-)linked to the :ref:`databases directory <ext_DBs>` (``db_path``)
subfolder ``filtering``, by default this is a FASTA file ``hg38.fa`` within the directory ``db_path/filtering``. 

IMP3 will remove rRNA reads from **metaT** reads. The user can decide which databases to use for this. The default is
to use the databases supplied by `SortMeRNA <https://bioinfo.lifl.fr/RNA/sortmerna/>`_. These databases should be available in 
the :ref:`databases directory <ext_DBs>` in FASTA format.

.. code-block:: yaml

  filtering:
    filter: hg38
  sortmerna:
    files:
      - rfam-5.8s-database-id98
      - silva-arc-16s-id95
      - silva-bac-16s-id90
      - silva-euk-18s-id95
      - rfam-5s-database-id98
      - silva-arc-23s-id98
      - silva-bac-23s-id98
      - silva-euk-28s-id98
	  
	  
--------
Assembly
--------

For the IMP3 :ref:`Assembly <step_assembly>` step, the user has the choice between two assemblers (``assembler``): `Megahit <http://www.metagenomics.wiki/tools/assembly/megahit>`_  
or `MetaSpades <http://cab.spbu.ru/software/meta-spades/>`_. Both assemblers are available for **metaG assemblies**. If the input are **metaG** and **metaT** reads, the user can 
choose between a **hybrid** assembly using both types of reads (``hybrid: true``, default) or a purely **metaG-only** assembly (``hybrid: false``). 
**Hybrid metaG/metaT assemblies** are currently only possible with `Megahit <http://www.metagenomics.wiki/tools/assembly/megahit>`_. `MetaSpades <http://cab.spbu.ru/software/meta-spades/>`_ is 
currently the only option to assemble **long and short metaG** reads together. 

For all assemblies, the user can set the minimal (``mink``) and maximal (``maxk``) *k*-mer size and the intervals between (``step``) for the multi 
*k*-mer assemblies performed by both assemblers. IMP3 performs :ref:`iterative assembly <step_assembly>` and the steps are merged using the 
`CAP3 <http://seq.cs.iastate.edu/cap3.html>`_ overlap assembler, which can be configured (``indentity``, ``overlap``).

All ``assembly`` config parameter are given below with default values:

.. code-block:: yaml
	  
  assembly: 
    hybrid: true
    assembler: megahit
    mink: 25
    maxk: 99
    step: 4
    cap3:
      identity: 98
      overlap: 100


-------------------------
Analysis: gene annotation
-------------------------

In the IMP3 :ref:`Analysis <step_analysis>` step, the predicted genes will be annotated using `hmmer <http://hmmer.org/>`_ and a set of HMM databases. 
The **HMM databases** must be present in the :ref:`databases directory <ext_DBs>`. The database name given in the ``hmm_DBs`` field must match the directory 
containing the HMMs.

For the IMP3 :ref:`Binning <step_binning>` step using :ref:`Binny <wtf_is_binny>`, contigs needs to have been annotated with the essential genes HMMs (``essential``).

Different HMM databases have different names and :ref:`pre-calibrated cut-offs <HMMs>`. If the HMM database has cut-off values from `hmmer <http://hmmer.org/>`_ 
calibration, set ``cutoff`` to ``--cut_tc``. Otherwise, keep an empty string.

Regarding the naming convention, some HMM databases contain multiple models for the same functional entity (e.g. multiple HMMs for one KEGG KO).
If the database contains models that are named with the functional entity + _ + another identifier, e.g. KO00033_1, the user can set
**trim** to ``--trimall`` and IMP3 will remove the underscore and last part of the identifier.

.. code-block:: yaml	  
	  
  hmm_DBs: "KEGG essential Pfam_A Resfams Cas dbCAN metacyc SwissProt TIGRPFAM"
  hmm_settings: 
    KEGG:
      cutoff: ""
      trim: "--trimall"
    essential:
      cutoff: "--cut_tc"
      trim: ""
    metacyc:
      cutoff: ""
      trim: "--trimall"
    Cas:
      cutoff: ""
      trim: ""
    Pfam_A:
      cutoff: "--cut_tc"
      trim: ""
    dbCAN:
      cutoff: ""
      trim: ""
    Resfams:
      cutoff: ""
      trim: ""
	
	
-----------------------------------------
Analysis: proteomics analysis preparation
-----------------------------------------

In the IMP3 :ref:`Analysis <step_analysis>` step, IMP3 provides several files for downstream **metaP** analysis.

- The user can use the protein sequence output of `prokka <http://www.vicbioinformatics.com/software.prokka.shtml>`_ directly as a search database for metaP analysis (not recommended). 
- IMP3 also provides a file with `Metaproteomics Analyzer <http://www.mpa.ovgu.de/>`_-conformant FASTA file headers, filtered to proteins having a minimum user-defined number of tryptic peptides (``filter_N_peptides``), with incomplete non-tryptic ends cut away. 
- IMP3 will add a FASTA file with host (or other reference) protein sequences (which should be placed in the :ref:`databases directory <ext_DBs>`). The provided host protein FASTA file should be consistent with the requirements of the proteomics software the user intends to use.
- IMP3 will optionally (``insert_variants: true/false``) generate a FASTA protein sequence file with the protein isoforms by inserting variants called on the assembled contigs (but :ref:`beware <variants_metaP>`).
	  
.. code-block:: yaml	  
	
  proteomics:
    filter_N_peptides: 2
    host_proteome: "hostproteomefile"
    insert_variants: false
  

------------------------------------------
Analysis: mapped reads per gene / function
------------------------------------------

In the IMP3 :ref:`Analysis <step_analysis>` step, `featureCounts <http://subread.sourceforge.net/>`_ is used count the numbers of **metaG** and/or 
**metaT** reads mapping to coding sequences and other annotated functional categories that can be used for downstream differential expression / abundance 
analyses. Per default, IMP3 assumes **metaG** reads to **not be strand-specific** ``mg: 0`` and **metaT** reads to be **stranded in the truSeq way** ``mt: 2``. 
	  
.. code-block:: yaml	  
	
  featureCountsStranding:
    mt: 2
    mg: 0  


--------
Taxonomy
--------

In the IMP3 :ref:`Taxonony <step_taxonomy> step, IMP3 runs `kraken2 <https://ccb.jhu.edu/software/kraken2/>`_ and `mOTUs2 <https://motu-tool.org/>`_ on the **metaG** reads. 
While the `mOTUs2 databases <https://motu-tool.org/>`_ are installed with the :ref:`IMP3 installation <install>`, the user needs to have a `kraken2 <https://ccb.jhu.edu/software/kraken2/>`_  
database available in the :ref:`databases directory <ext_DBs>` and give its name under ``krakendb``.

If the :ref:`Taxonony <step_taxonomy>` and :ref:`Analysis <step_analysis>` steps are performed, IMP3 will search 
the predicted genes by default for the mOTUs2-marker genes. The default COGs can be customized (field ``COGS``).
	  
.. code-block:: yaml

  krakendb: minikraken2
  COGS: "COG0012 COG0018 COG0215 COG0525 COG0541 COG0016 COG0172 COG0495 COG0533 COG0552"

-------
Binning
-------

In the IMP3 :ref:`Binning <step_binning>` step, the user can choose to run any combination of three different ``binners``: `MaxBin2 <https://sourceforge.net/projects/maxbin2/>`_, 
`MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>`_, :ref:`Binny <wtf_is_binny>`.
The results will then refined using `DAStool <https://github.com/cmks/DAS_Tool>`_.

For `MaxBin2 <https://sourceforge.net/projects/maxbin2/>`_ and :ref:`Binny <wtf_is_binny>`, **metaG** reads are required.
If there are only **metaT** reads, only `MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>`_ will run.

Each binner comes with a few parameters that the user can configure:

- The minimum contig lengths ``cutoff`` for `MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>`_ is **1500**. 
- The default minimum contig lengths ``cutoff`` for `MaxBin2 <https://sourceforge.net/projects/maxbin2/>`_ is **1000**.
- :ref:`Binny <wtf_is_binny>` uses the minimum contig lengths cutoff of ``vizbin`` (default **1000**).
- If :ref:`VizBin <http://claczny.github.io/VizBin/> output should be included into the summary, :ref:`Binny <wtf_is_binny>` has to be chosen.

.. code-block:: yaml

  binning:
    binners: "MaxBin MetaBAT binny"
    MaxBin:
      cutoff: 1000
    MetaBAT:
      cutoff: 1500
    binny:
      pk: 10
      nn: 4
      cutoff: 1000
    vizbin:
      dimension: 50
      kmer: 5
      size: 4
      theta: 0.5
      perp: 30
      cutoff: 1000

