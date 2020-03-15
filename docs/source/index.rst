.. IMP3 documentation master file, created by
   sphinx-quickstart on Sat Mar 14 17:06:18 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

IMP3
====

The Integrated Meta-Omics Pipeline, IMP, is a reproducible and modular pipeline for large-scale 
integrated analysis of coupled metagenomic and metatranscriptomic data. IMP incorporates read 
preprocessing, assembly (including optional iterative co-assembly of metagenomic and 
metatranscriptomic data), analysis of microbial community structure (taxonomy) and function 
(gene functions and abundances), binning, as well as summaries and visualizations. 
As an added functionality, IMP also performs single-omic metagenomic or metatranscriptomic analyses.

IMP is implemented in `Snakemake <https://snakemake.readthedocs.io/en/stable/index.html>`_.


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
