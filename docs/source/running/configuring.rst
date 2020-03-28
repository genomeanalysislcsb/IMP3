.. _configuration:

=====================
Configuring IMP3 runs
=====================

You can choose the steps IMP3 should perform, the data it should run on, and the settings
of each step by means of a `.yaml <https://yaml.org/>`_ file. You have to point snakemake to the config file when you :ref:`run <run_IMP>`
IMP3, using the ``--configfile`` argument.

This description will guide you through all sections of the config file, but note that the order of the items in the config file is arbitrary.
You can also omit most fields and IMP3 will use default settings. Exceptions are the :ref:`inputs <input_options>` and
an output directory. IMP3 will write the used configuration including the defaults into the output folder as ``sample.config.yaml``.

----------------------
Inputs and directories
----------------------

You need to tell IMP3 where the **raw or input files** are located. (The default is no input, and nothing will be done). Give absolute paths. 
You also need to tell IMP3 where the **ouput** should be written. Please also give an absolute path here.
It's advisable to give a **sample name** to IMP3. Avoid special characters. Punctuation will be converted to underscore.
The sample name will be prepended to, e.g. the contig names. (If you don't give a sample name, IMP3 will take the last two parts of the output
path, concatenate them with an underscore and use the results as sample name).
A directory for **temporary files** can be given. If you don't give an absolute path, it will be put into the output directory
(and, for the default, named tmp).
If you let IMP3 perform any steps that need :ref:`databases <ext_DBs>` (i.e. preprocessing (filtering of reads against a host genome, filtering rRNA reads),
analysis (gene calling, functional annotation), or taxonomy (kraken2)), you also need to give the path to one directory that houses all the
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
a space under **steps**. The default is to run all steps.
If you choose to perform pre-processing of reads, you can turn on or off **preprocessing_filtering** (i.e. the removal of reads which
map against a reference genome, most commonly a host genome). The default is for filtering to be on and reads to be mapped against
a file hg38.fa in the database directory.
If you choose to let IMP3 summarize your results, you can define which steps should be taken. The default is to extract some
**stats** on all performed steps and to visualize them (**vis**). You can choose to do neither, in which case IMP3 will not make summaries.
You can have only stats, or stats and visualization by stating one or both separated by a space.
Don't expect IMP3 to do anything useful, if you choose only vis (visualization). 

.. code-block:: yaml

  steps: "preprocessing assembly analysis binning taxonomy summary"
  preprocessing_filtering: true
  summary_steps: "stats vis"

.. _hardware:

------------------
Available hardware 
------------------

Tell IMP3 how large the **RAM** of your **normal** and **bigmem** cpus is and how many **cores** of each you want to maximally use in a step.
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


