
.. _output_analysis:

=====================
IMP3 output: Analysis
=====================

----------
Annotation
----------

Final annotations:	annotation_CDS_RNA_hmms.gff
Feature counts funct.:	<ome>.<database>.counts.<tsv|tsv.summary>
Feature counts genes:	<ome>.<CDS|rRNA>_counts.<tsv|tsv.summary>
Proteomics DB:	proteomics.final.faa
intermediary:	 proteomics.proteins.faa
Prokka results:	annotation.filt.gff
sequences:	prokka.<faa|ffn|fna|fsa>
indices:	prokka.<faa|ffn|fna|fsa>.<suffix>
others:	prokka.<suffix>
HMMer intermediary:	intermediary.tar.gz
DASTool-ID:	 annotation_CDS_RNA_hmms. contig2ID.tsv

SNPs:		<ome>.variants.samtools.vcf.gz
Indices:	<ome>.variants.samtools.vcf.gz.tbi

----
SNPs
----

--------
Taxonomy
--------

The **Analysis** directory contains the output of the :ref:`taxonomy <output_taxonomy>` step. The exception
is the output of GTDBtk, which is run on the bins selected by DASTool. See the :ref:`next <output_taxonomy>`
section for details.

.. _output_analysis_stats:

----------------------------
Stats from the analysis step
----------------------------

In addition to the annotation and SNP calling, the analysis step will also collect some stats for you and save them to the ``Stats``
directory. 

