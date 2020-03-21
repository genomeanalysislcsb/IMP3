.. _configuration:

=====================
Configuring IMP3 runs
=====================

You can choose the steps IMP3 should perform, the data it should run on, and the settings
of each step by means of a .`yaml <https://yaml.org/>`_ file. 

This description will guide you through all sections, but note that the order of the items in the config file is arbitrary.
You can also omit most fields and IMP3 will use default settings. Exceptions are the :ref:`inputs <input_options>` and
an output directory.


----------------------
Inputs and directories
----------------------

You need to tell IMP3 where the raw or input files are located. Give absolute paths. (The default is no input, and nothing will be done.)
You also need to tell IMP3 where the ouput should be written. Please also give an absolute path here.
It's advisable to give a sample name to IMP3. This will be prepended to, e.g. the contig names.
A directory for temporary files can be given. If you don't give an absolute path, it will be put into the output directory (and named tmp).
If you let IMP3 perform any steps that need databases (i.e. preprocessing (filtering of reads against a host genome, filtering rRNA reads),
analysis (gene calling, functional annotation), or taxonomy (kraken)), you also need to give the path to one directory that houses all the
databases.

.. code-block:: yaml

  raws: 
    Metagenomics: ""
    Metatranscriptomics: ""
    LongReads: ""
    LongReadTech: ""
    Contigs: ""
    Alignment_metagenomics: ""
    Alignment_metatranscriptomics: ""
  outputdir: ""
  sample: test
  tmp_dir: tmp
  db_path: "/data/idiv_biu/heintzbu/EXTDB/IMP_DBs"


----------
IMP3 steps
----------

With IMP3, you can choose to run the whole workflow or just some :ref:`steps <steps_overview>`. List the steps, separated by
a space under ``steps``. The default is to run all steps.
If you choose to perform pre-processing of reads, you can turn on or off filtering (i.e. the removal of reads which
map against a reference genome, most commpnly a host genome). The default is for filtering to be on and reads to be mapped against
a file hg38.fa in the database directory.
If you choose to let IMP3 summarize your results, you can define which steps should be taken. The default is to extract some
stats on all performed steps and to visualize them. You can choose to do neither, in which case IMP3 will not make summaries.
You can have only stats, or stats and visualization by stating one or both separate by a space.
Don't expect IMP3 to do anything useful, if you choose only vis (visualization). 

.. code-block:: yaml

  steps: "preprocessing assembly analysis binning taxonomy summary"
  preprocessing_filtering: true
  summary_steps: "stats vis"

------------------
Available hardware 
------------------

Tell IMP3 how large the RAM of you normal and bigmem cpus is and how many cores you want to maximally use in a step.

.. code-block:: yaml

  mem:
    normal_mem_per_core_gb: 4
    big_mem_total_gb: 1600
    big_mem_cores: 8
    big_mem_per_core_gb: 200







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
assembly: 
  hybrid: true
  assembler: megahit
  mink: 25
  maxk: 99
  step: 4
  cap3:
    identity: 98
    overlap: 100
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
  SwissProt:
    cutoff: ""
    trim: "--trimall"
  TIGRPFAM:
    cutoff: ""
    trim: ""
  dbCAN:
    cutoff: ""
    trim: ""
  Resfams:
    cutoff: ""
    trim: ""
COGS: "COG0012 COG0018 COG0215 COG0525 COG0541 COG0016 COG0172 COG0495 COG0533 COG0552"
prokka: 
  prefix: prokka
  databases: 
    - "cm/Bacteria.i1i"
    - "genus/Staphylococcus.phr"
    - "hmm/CLUSTERS.hmm.h3f"
    - "kingdom/Archaea/sprot.phr"
featureCountsStranding:
  mt: 2
  mg: 0
proteomics:
  filter_N_peptides: 2
  host_proteome: ""
  insert_variants: false
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
krakendb: minikraken2

