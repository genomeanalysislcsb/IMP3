.. _configuration:

=====================
Configuring IMP3 runs
=====================

IMP3 is designed modular and highly configurable. Users can define the analysis steps IMP3 should perform, the data it should run on, the databases and tools 
that should be used and the settings of each step in a single **config file** in `yaml <https://yaml.org/>`_ format. This **config file** is provided `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ 
as described here :ref:`run <run_IMP>`, using the ``--configfile`` argument.

This description will guide you through all sections of the **config file**. 

Some general notes on the **config file**:

- the order of the items in the config file is arbitrary.
- the **config file** does not need to contain all fields.
- only two fields are required: **:ref:`inputs <input_options>`** and an **output directory**. 
- all other fields are optional and if emtpy the IMP3 default settings will be used. 
- IMP3 will write the user-specified configuration including pre-defined fields the output folder as ``sample.config.yaml``.

-------------------------------
Inputs, outputs and directories
-------------------------------

As required input the user has to tell IMP3 where the **raw or input files** are located. (The default is no input, and nothing will be done). 
Different steps in IMP3 require different input files. All input files are given with their absolute paths. Input files can be:

- **MetaG** and/or **MetaT** paired-end or single-end raw `Fastq <https://en.wikipedia.org/wiki/FASTQ_format>`_ files.
- **MetaG** and/or **MetaT** `BAM <https://genome.sph.umich.edu/wiki/BAM>`_ files (requires ``Contigs`` file, required if the :ref:`Preprocessing <step_preprocessing>` and :ref:`Assembly <step_assembly>` steps are skipped).
- ``Contigs`` file in FASTA format (used for the mapping of the metaG und metaG raw files, required if the :ref:`Preprocessing <step_preprocessing>` and :ref:`Assembly <step_assembly>` steps are skipped).
- ``LongReads`` file in `Fastq <https://en.wikipedia.org/wiki/FASTQ_format>`_ (requires a ``LongReadTech`` description).
- ``LongReadTech`` description of used long read sequencing technology.
- ``Gff`` file to gff based on the ``Contigs`` file (only required if the :ref:`Preprocessing <step_preprocessing>`, :ref:`Assembly <step_assembly>` and :ref:`Analysis <step_analysis>` steps are skipped and only the step :ref:`Binning <step_binning>` should performed).

Paired end read files are separated by ``<space>``.

The ``outputdir`` (with full path) defines the directory, where the IMP3 workflow output should be written to.

Users can give a **sample name**: ``sample`` (default value is ``test``). Special characters should be avoided, punctuation will be converted to underscore.
The **sample name** will be used throughout the IMP3 workflow, e.g. prepended to contig names. If the user-defined sample name is empty (``sample: ""``), IMP3 
will take the last two parts of the output path, concatenate them with an underscore and use the results as sample name).

A directory for **temporary files** can be given (by default it is named ``tmp``). If you don't give an absolute path, it will be put into the output 
directory.

If you let IMP3 perform any steps that need :ref:`databases <ext_DBs>` (i.e. :ref:`Preprocessing <step_preprocessing>` (filtering of reads against a host genome, filtering rRNA reads),
:ref:`Analysis <step_analysis>` (gene calling, functional annotation), or :ref:`Taxonomy <step_taxonomy>` (kraken2)), you also need to give the full path to a single directory that houses all the
databases.

.. code-block:: yaml

  raws: 
    Metagenomics: "/path/to/pair1 /path/to/pair2"
    Metatranscriptomics: "/path/to/pair1 /path/to/pair2"
    LongReads: "/path/to/longreads"
    LongReadTech: ""
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

With IMP3, you can choose to run the whole workflow or just some :ref:`steps <steps_overview>`. The field ``steps`` lists the user-defined IMP3 analysis steps, separated by
a space, under **steps**. The default is to run all steps.

By default (``preprocessing_filter: true``) and having raw reads as input, IMP3 will pre-process the raw sequencing reads (i.e. the removal of reads 
which map against a reference genome, most commonly a host genome). The user can set ``preprocessing_filtering`` to ``false`` to skip this step. 
The default filtering will quality-trim reads and map them mapped against the human genome given as file ``hg38.fa`` in the database directory under ``db_path``.

If ``steps`` contains ``summary``, IMP3 will summarize your results. The user can define which steps should be taken for summary in the field ``summary_steps``. 
The default is to extract summary statistics (``stats``) on all performed steps and to visualize them (``vis``). If ``summary_steps`` is empty, IMP3 will 
not make summaries. Currently, the user can either choose ``stats`` or ``vis`` or both. 

.. code-block:: yaml

  steps: "preprocessing assembly analysis binning taxonomy summary"
  preprocessing_filtering: true
  summary_steps: "stats vis"

.. _hardware:

--------------------------
Available hardware support
--------------------------

Some IMP3 steps need high memory machines depending on the size of the dataset and the databases used. IMP3 is able to run tasks on specific 
type of nodes, depending on the **memory** demands of the jobs, high memory demanding jobs are automatically scheduled to ``big_mem`` nodes. 
Depending on how much **RAM** your local computer or compute cluster have available for your jobs the user can define the settings for  **normal** 
and **bigmem** cpus:

- ``normal_mem_per_core_gb`` (max. available RAM (in Gb) per core of **normal** compute cores, default 4 Gb)
- ``big_mem_total_gb`` (max. available total RAM (in Gb) of **bigmem** compute cores, default 1600 Gb)
- ``big_mem_cores`` (available number of **bigmem** compute cores, default 8)
- ``big_mem_per_core`` (max. available RAM (in Gb) per core of **bigmem** compute cores, default 200 Gb)

The number of cores is determined by the ``--cores`` argument of your `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ :ref:,call <run_IMP>`, rather than set here.

The jobs will be scheduled depending on the availability and RAM demands depending on how many **cores** of each the user maximally use in a step.

*Defaults may be blatantly inappropriate for your system. Note that if you don't let IMP3 submit jobs for you, most steps are run independently of how much RAM is available, but some do (e.g. mapping, kraken).*

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

For the IMP3 :ref:`Preprocessing <step_preprocessing>` step, IMP3 is using **:ref:`Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`** for read trimming.
The user can set all **:ref:`Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`** options (defaults given below):

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

If your reads originate from an :ref:`Illumina Nextseq <https://www.illumina.com/systems/sequencing-platforms/nextseq.html?langsel=/us/>` machine, we advise to set **nextseq** to true. This will remove
trailing Gs that the Nextseq systems add to reads that are shorter than the run cycles.

-------------------------
Pre-processing: filtering 
-------------------------

The IMP3 :ref:`Preprocessing <step_preprocessing>` step allows to map against a user-defined reference genome (e.g. host). 
The user can set the file name of the reference genome (``filter``, which should be located or linked to the :ref:`databases directory <ext_DBs>` (``db_path``)
subfolder ``filtering``, e.g. this would be a Fasta file ``hg39.fa`` within the directory ``db_path/filtering``. 

IMP3 will remove rRNA reads from your **metaT** reads. The user can decide which databases to use for this. The default is
to use the databases supplied by :ref:`SortMeRNA <https://bioinfo.lifl.fr/RNA/sortmerna/>`. These databases should be available in 
the :ref:`databases directory <ext_DBs>` in Fasta format.

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

For the IMP3 :ref:`Assembly <step_assembly>` step, the user has the choice between two assemblers, :ref:`megahit <http://www.metagenomics.wiki/tools/assembly/megahit>` 
or :ref:`Metaspades <http://cab.spbu.ru/software/meta-spades/>` as **metaG assemblers**. If you have **metaG** and **metaT** reads, the user can 
choose between a **hybrid** assembly using both types of reads or a purely **metaG-only** assembly. 

**Hybrid metaG/metaT assemblies** are currently only possible with :ref:`megahit <http://www.metagenomics.wiki/tools/assembly/megahit>`. 

:ref:`Metaspades <http://cab.spbu.ru/software/meta-spades/>` is currently your only option, if you want to assemble **long and short metaG** reads together.

For all assemblies, the user can set the minimal and maximal *k*-mer size and the intervals between for the multi *k*-mer assemblies performed by both assemblers. 
IMP3 performs :ref:`iterative assembly <step_assembly>` and the steps are merged using the :ref:`Cap3 <http://seq.cs.iastate.edu/cap3.html>` overlap assembler, which you can configure.
	  
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


----------------------
Analysis: gene calling
----------------------

In the IMP3 :ref:`Analysis <step_analysis>` step, IMP3 perfoms gene calling using a customized version of :ref:`Prokka <http://www.vicbioinformatics.com/software.prokka.shtml>` 
on the assembled contigs (see :ref:`customized prokka <prokkaC>`) to call open reading frames, rRNA and tRNA genes. The default databases of prokka 
should be placed in your :ref:`databases directory <ext_DBs>` and are part of the default configuration.
	  
.. code-block:: yaml

  prokka: 
    databases: 
      - "cm/Bacteria.i1i"
      - "genus/Staphylococcus.phr"
      - "hmm/CLUSTERS.hmm.h3f"
      - "kingdom/Archaea/sprot.phr"
	  
	  
-------------------------
Analysis: gene annotation
-------------------------

In the IMP3 :ref:`Analysis <step_analysis>` step, the called open reading frames will be annotated using :ref:`hmmer <http://hmmer.org/>` and a set of HMM databases. 
The **HMM databases** must be present in your :ref:`databases directory <ext_DBs>`. The name given in the ``hmm_DBs`` field must match the directory 
containing the HMMs.

For the IMP3 :ref:`Binning <step_binning>` step using :ref:`Binny <wtf_is_binny>`, contigs needs to be annotated with the essential genes HMMs (``essential``).

Different HMM databases have different names and :ref:`pre-calibrated cut-offs <HMMs>`. If your HMM database has cut-off values from :ref:`hmmer <http://hmmer.org/>` 
calibration, set ``cutoff`` to ``--cut_tc``. Otherwise, keep an empty string.

Regarding the naming, some HMM databases contain multiple models for the same functional entity (e.g. multiple HMMs for one KEGG KO).
If you happen to have such a collection, and the models are named with the functional entity + _ + another identifier, you can set
**trim** to ``--trimall`` and IMP3 will remove the underscore and final identifier for you.
	  
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
	
	
--------------------------------
Analysis: proteomics preparation
--------------------------------

In the IMP3 :ref:`Analysis <step_analysis>` step, IMP3 provides several files for downstream **metaP** analysis.

- The user can the protein sequence output of :ref:`Prokka <http://www.vicbioinformatics.com/software.prokka.shtml>` as a search database for metaP analysis. 
- IMP3 also provides a file with :ref:`Metaproteomics Analyzer <http://www.mpa.ovgu.de/>`-conformant Fasta file headers, filtered to proteins having a 
minimum user-defined number of tryptic peptides (``filter_N_peptides``), with incomplete non-tryptic ends cut away. 
- IMP3 will add a Fasta file with host (or other reference) protein sequences (which should be placed in your :ref:`databases directory <ext_DBs>`). 
- IMP3 will generates optional (``insert_variants: true/false``) a Fasta protein file with the protein isoforms by inserting variants called on the assembled contigs (see :ref:`beware <variants_metaP>`).
	  
.. code-block:: yaml	  
	
  proteomics:
    filter_N_peptides: 2
    host_proteome: "hostproteomefile"
    insert_variants: false
  

-----------------------------------
Analysis: reads per function / gene
-----------------------------------

In the IMP3 :ref:`Analysis <step_analysis>` step, IMP3 runs :ref:`featureCounts <http://subread.sourceforge.net/>` to count the numbers of **metaG** and/or 
**metaT** reads mapping to coding sequences and other annotated functional categories that can be used for downstream differential expression / abundance 
analyses. Per default, IMP3 assumes **metaG** reads to **not be strand-specific (mt: 0)** and **metaT** reads to be **stranded in the truSeq way (mt: 2)**. 
	  
.. code-block:: yaml	  
	
  featureCountsStranding:
    mt: 2
    mg: 0  


--------
Taxonomy
--------

In the IMP3 :ref:`Taxonony <step_taxonomy>` step, IMP3 runs :ref:`kraken2 <https://ccb.jhu.edu/software/kraken2/>` and :ref:`mOTUs2 <https://motu-tool.org/>` on the **metaG** reads. 
While the :ref:`mOTUs2 databases <https://motu-tool.org/>` are installed with the :ref:`IMP3 installation <install>`, the user need to have a :ref:`kraken2 <https://ccb.jhu.edu/software/kraken2/>` 
database available in your :ref:`databases directory <ext_DBs>` and give its name under ``krakendb``.

If the IMP3 input are **metaG** reads and the :ref:`Taxonony <step_taxonomy>` and :ref:`Analysis <step_analysis>` steps are performed, IMP3 will search 
the predicted genes for for the mOTUs2-marker genes. The default COGs can be customized (field ``COGS``).
	  
.. code-block:: yaml

  krakendb: minikraken2
  COGS: "COG0012 COG0018 COG0215 COG0525 COG0541 COG0016 COG0172 COG0495 COG0533 COG0552"


-------
Binning
-------

In the IMP3 :ref:`Binning <step_binning>` step, the user can choose to run any combination of three different ``binners``: :ref:`MaxBin2 <https://sourceforge.net/projects/maxbin2/>`, 
:ref:`MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>`, :ref:`Binny <wtf_is_binny>`.
The results of the all binnings will then refined using :ref:`DAStool <https://github.com/cmks/DAS_Tool>`.

For :ref:`MaxBin2 <https://sourceforge.net/projects/maxbin2/>` and :ref:`Binny <wtf_is_binny>`, it is required that you have **metaG** reads.
If there are only **metaT** reads / assemblies are available, only :ref:`MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>` will run.

Each binner comes with a few parameters that the user can configure:

- The default minimum contig lengths ``cutoff`` for :ref:`MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/src/master/>` is *1500*. 
- The default minimum contig lengths ``cutoff`` for :ref:`MaxBin2 <https://sourceforge.net/projects/maxbin2/>` is *1000*.
- :ref:`Binny <wtf_is_binny>` has no minimum ontig lengths cutoff. If :ref:`VizBin <http://claczny.github.io/VizBin/>` output should be included into the summary, 
:ref:`Binny <wtf_is_binny>` has to be chosen.

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


