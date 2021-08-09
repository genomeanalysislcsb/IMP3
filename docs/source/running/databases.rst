.. _ext_DBs:

==================
External databases
==================

IMP3 makes use of external databases. They should be provided in or linked to one directory. See the :ref:`install <installation guide>`.

The following steps access external databases:

* Preprocessing: 

  * adapter file for trimmomatic, 

  * one or more genomes to filter against (optional), 

  * sortmeRNA Silva database (metatranscriptomics only).

* Taxonomy:

  * a kraken2 database,

  * the GTDB-tk database (only, if binning is done),

  * the eukdetect (BUSCO) database (optional).

* Annotation:

  * HMM collections.

