
.. _step_analysis:

====================
IMP3 steps: Analysis
====================

In the Analysis step, IMP3 calls open reading frames, rRNA and tRNA genes, and annotate CRISPR repeats.
Open reading frames are functionally annotated. The number of reads mapping to each gene and functional
group of genes are calculated. The steps are described in more detail below.

.. _prokkaC:

-----------------
Customized prokka
-----------------

IMP3 uses `prokka <http://www.vicbioinformatics.com/software.prokka.shtml>`_ to call 
open reading frames (ORF), rRNA and tRNA genes, and annotate CRISPR repeats. Internally,
`prokka <http://www.vicbioinformatics.com/software.prokka.shtml>`_ calls `prodigal <https://github.com/hyattpd/Prodigal>`_
for the ORF calling, `barrnap <https://github.com/tseemann/barrnap>`_ for the rRNA regions,
`ARAGORN <https://academic.oup.com/nar/article/32/1/11/1194008>`_ for tRNA loci and
`MinCED <https://github.com/ctSkennerton/minced>`_ to detect CRISPR arrays.

Prokka forces prodigal to only call complete genes. Due to the fragmented nature of metagenomic contigs, it
is preferable to also allow partial genes. IMP3's customized prokka allows prodigal to call incomplete 
ORFs and records whether prodigal detected start and stop codons. One side-effect of this is that the
amino acid sequences prokka returns ``prokka.faa`` are badly formatted. The prokka amino acid sequences also don't start with M
if prodigal called genes with an alternative start codon. Both issues are corrected in another IMP3 step as part of the metaproteomics
preparations (``proteomics.proteins.faa``).

Prokka would usually provide some functional analyses by aligning the called ORFs to some databases. However, this analysis
is optimized for speed, meaning that genes that have been annotated with one database are not annotated with the next. This
leads to genes potentially having inconsistent annotations, and it would be impossible for the user to find out what would have
been the best hit in another database. Since IMP3 does :ref:`functional annotations <HMMs>` of all genes with any database
the user chooses to reach consistent annotations, we've disabled the prokka-based annotation.

Prokka also spends considerable time to convert its output into genbank format. As IMP3 has no need for genbank-formatted
data, we've disabled this.

.. _HMMs:

--------------------
Annotation with HMMs
--------------------





.. _variants_metaP:

---------------------------
Variants for metaP searches
---------------------------
