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

IMP is implemented in _Snakemake: https://snakemake.readthedocs.io/en/stable/index.html .


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

   running/configuring
   running/run





.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
