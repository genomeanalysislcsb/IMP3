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
The **sample name** will be used throughout the IMP3 workflow, e.g. prepended to contig names. If the user-defined sample name is empty (``sample=""``), IMP3 
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

IMP3 is able to run tasks on specific type of nodes, depending if the user d how large the **RAM** of your **normal** and **bigmem** cpus is and how many **cores** of each you want to maximally use in a step.
*Defaults may be blatantly inappropriate for your system. Note that if you don't let IMP3 submit jobs for you,
most steps don't care about how much RAM you say you have, but some do (e.g. mapping, kraken).*
The number of normal cores is determined by the ``--cores`` argument of your snakemake :ref:,call <run_IMP>`, rather than set here.
Currently you need to have both kinds of cores, if you run IMP3 in :ref:`cluster mode <run_cluster>`,
or you have to change the :ref:`cluster config file <cluster_config>`. 

.. code-block:: yaml

  mem:
    normal_mem_per_core_gb: 4
    big_mem_total_gb: 1600
    big_mem_cores: 8
    big_mem_per_core_gb: 200


------------------------
Pre-processing: trimming 
------------------------

If you want IMP3 to pre-process your reads, you can set all **trimmomatic** options (defaults below).
If your reads originate from an Illumina Nextseq machine, we advise to set **nextseq** to true. This will remove
trailing Gs that the Nextseq systems add to reads that are shorter than the run cycles.

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

If you want IMP3 to pre-process your reads and you've chosen to map against a reference genome (e.g. host),
you can set the file name of the reference genome, which should be located or linked to in the :ref:`databases directory's <ext_DBs>`
subfolder **filtering**.
IMP3 will remove rRNA reads from your metatranscriptomic reads. You can decide which databases to use for this. The default is
to use the databases supplied by sort me RNA. You should have these in the :ref:`databases directory <ext_DBs>`, too.

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

If you let IMP3 do an **assembly**, you have the choice between Megahit and Metaspades as metagenome **assemblers**.
If you have metagenomic and metatranscriptomic reads, you can choose between a **hybrid** assembly of both or a purely metagenomic
assembly. Hybrid metaG/metaT assemblies are currently only possible with Megahit. On the other hand,
Metaspades is your only option, if you want to assemble long and short metagenomic reads.
You can set the minimal and maximal *k*-mer size and the intervals between. IMP3 performs an :ref:`iterative assembly <step_assembly>`
and the steps are merged using the Cap3 overlap assembler, which you can configure.
	  
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

If you want IMP3 to analyse the assembled contigs, it will run a :ref:`customized prokka <prokkaC>` to call open reading frames,
rRNA and tRNA genes.
The default databases of prokka should be in your :ref:`databases directory <ext_DBs>` and are part of the default configuration.
	  
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

The called open reading frames will be attempted to be annotated using HMMs. The HMMs must be present in your :ref:`databases directory <ext_DBs>`.
The name given in the **hmm_DBs** field must match the directory containing the HMMs.
If you want IMP3 to later perform binning using Binny, you should let IMP3 at least annotated the essential genes.
Different HMM collections have different names and :ref:`pre-calibrated cut-offs <HMMs>`. If your collection has cut-off values from calibration,
set **cutoff** to ``--cut_tc``. Otherwise, keep an empty string.
Regarding the naming, some HMM collections contain multiple models for the same functional entity (e.g. multiple HMMs for one KEGG KO).
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

You could just use the amino acid output of prokka as a database for metaproteomics. IMP3 also gives you a file with
Metaproteomics Analyzer-conformant FASTA headers, filtered to proteins having at least as many tryptic peptides as you
wish, with incomplete non-tryptic ends cut away. It will also add a file with host (or other reference) proteins (which should be
present in your :ref:`databases directory <ext_DBs>`), if you give this file's name. Finally, it can insert variants, but :ref:`beware <variants_metaP>`.
	  
.. code-block:: yaml	  
	
  proteomics:
    filter_N_peptides: 2
    host_proteome: ""
    insert_variants: false
  

-----------------------------------
Analysis: reads per function / gene
-----------------------------------

IMP3 runs featureCounts to give you numbers of reads mapping to coding sequences and functional categories. You may want to use these
for differential expression / abundance analyses. Per default, IMP3 assumes metagenomics reads to **not be strand-specific (0)** and
metatranscriptomics reads to be **stranded in the truSeq way (2)**. 
	  
.. code-block:: yaml	  
	
  featureCountsStranding:
    mt: 2
    mg: 0  


--------
Taxonomy
--------

If you give reads to IMP3 and choose taxonomy to be done, IMP3 will run kraken2 and mOTUs2 on your reads. While the mOTUs databases
are installed with the :ref:`IMP3 install <install>`, you need to have a kraken2 database in your :ref:`databases directory <ext_DBs>`
and give its name under **krakendb**.
If you give reads to IMP3 and ask for taxonomy and analysis to be done, it will also go looking for the mOTU-marker genes in your
predicted genes and let you know which ones it found. The default COGs can be adjusted, if you feel adventurous.
	  
.. code-block:: yaml

  krakendb: minikraken2
  COGS: "COG0012 COG0018 COG0215 COG0525 COG0541 COG0016 COG0172 COG0495 COG0533 COG0552"


-------
Binning
-------

If you ask IMP3 to bin your contigs, you can run MaxBin2 and/or MetaBAT2 and/or :ref:`binny <wtf_is_binny>` on them, prior to DAStool.
For MaxBin2 and binny, it is required that you have metagenomics depths: if you have only metatranscriptomics reads / assemblies,
only MetaBAT2 will run.
Each binner comes with a few settings that you may want to configure.
The minimum cutoff for MetaBAT is 1500. binny has no minimum number.
If you want to get VizBin output, you have to choose :ref:`binny <wtf_is_binny>`.

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


