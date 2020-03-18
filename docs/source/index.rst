IMP3
====

The Integrated Meta-Omics Pipeline, IMP, is a reproducible and modular pipeline for large-scale 
integrated analysis of coupled metagenomic and metatranscriptomic data. IMP incorporates read 
preprocessing, assembly (including optional iterative co-assembly of metagenomic and 
metatranscriptomic data), analysis of microbial community structure (taxonomy) and function 
(gene functions and abundances), binning, as well as summaries and visualizations. 

IMP is implemented in `Snakemake <https://snakemake.readthedocs.io/en/stable/index.html>`_. It is user-friendly
and can be adjusted easily configured to perform single-omic metagenomic or metatranscriptomic analyses, as
well as single modules, e.g. only assembly or only binning.

IMP is being developed at the `Luxembourg Centre for Systems Biomedicine <https://wwwen.uni.lu/lcsb/research/eco_systems_biology>`_ and the
`German Centre for Integrative Biodiversity Research (iDiv) Halle-Jena-Leipzig <https://www.idiv.de/en/groups_and_people/central_management/bioinformatics_unit_biu/research.html>`_.


.. toctree::
   :maxdepth: 1
   :caption: Installation
   :hidden:
   :name: installation

   installation/installation

.. toctree::
   :maxdepth: 1
   :caption: Configuring and running IMP
   :hidden:
   :name: running_IMP

   running/run
   running/configuring
   running/inputs
   running/cluster_config
   
.. toctree::
   :maxdepth: 1
   :caption: IMP's steps
   :hidden:
   :name: steps

   steps/Overview
   steps/Preprocessing
   steps/Assembly
   steps/Analysis
   steps/Binning
   steps/Taxonomy
   steps/Summary
   
.. toctree::
   :maxdepth: 1
   :caption: IMP's output
   :hidden:
   :name: outputs

   output/Overview
   output/Report
   output/Stats
   output/Preprocessing
   output/Assembly
   output/Analysis
   output/Binning
   output/Taxonomy
   
.. toctree::
   :maxdepth: 1
   :caption: Configuring and running multi-IMP
   :hidden:
   :name: running_multiIMP

   multi_running/run
   multi_running/configuring
   multi_running/sample_file
   
.. toctree::
   :maxdepth: 1
   :caption: multi-IMP's steps
   :hidden:
   :name: multi_steps

   multi_steps/Overview
   multi_steps/Annotation
   multi_steps/Taxonomy
   multi_steps/Catalogue
   multi_steps/Dereplication
   multi_steps/Re-binning
   multi_steps/Summary
   
.. toctree::
   :maxdepth: 1
   :caption: multi-IMP's output
   :hidden:
   :name: multi_outputs

   multi_output/Overview
   multi_output/Report
   multi_output/Stats
   multi_output/Annotation
   multi_output/Taxonomy
   multi_output/Catalogue
   multi_output/Dereplication
   multi_output/Re-binning


.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
