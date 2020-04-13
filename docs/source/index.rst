####
IMP3
####

The Integrated Meta-Omics Pipeline, IMP, is a reproducible and modular pipeline for large-scale 
integrated analysis of coupled metagenomic and metatranscriptomic data. IMP3 incorporates read 
preprocessing, assembly (including optional iterative co-assembly of metagenomic and 
metatranscriptomic data), analysis of microbial community structure (taxonomy) and function 
(gene functions and abundances), binning, as well as summaries and visualizations. 

IMP3 is implemented in `Snakemake <https://snakemake.readthedocs.io/en/stable/index.html>`_. It is user-friendly
and can be easily configured to perform single-omic metagenomic or metatranscriptomic analyses, as
well as single modules, e.g. only assembly or only binning.

IMP3 is being developed at the `Luxembourg Centre for Systems Biomedicine <https://wwwen.uni.lu/lcsb/research/eco_systems_biology>`_ and the
`German Centre for Integrative Biodiversity Research (iDiv) Halle-Jena-Leipzig <https://www.idiv.de/en/groups_and_people/central_management/bioinformatics_unit_biu/research.html>`_.


.. toctree::
   :maxdepth: 1
   :caption: Installation
   :hidden:
   :name: installation

   installation/installation
   
.. toctree::
   :maxdepth: 1
   :caption: IMP3
   :hidden:
   :name: IMP3
   
   running/run_idx
   steps/Overview
   output/Overview
   
.. toctree::
   :maxdepth: 1
   :caption: Multi-IMP3
   :hidden:
   :name: running_multiIMP3

   multi_running/run_idx
   multi_steps/Overview
   multi_output/Overview

.. toctree::
   :maxdepth: 1
   :caption: FAQ
   :hidden:
   :name: FAQ

   FAQ

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
