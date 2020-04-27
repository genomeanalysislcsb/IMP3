
.. _step_preprocessing:

=========================
IMP3 steps: Preprocessing
=========================

Preprocessing of reads consists of 1 to 3 steps: trimming,
removal of ribosomal RNAs, filtering of reads mapping to one or more reference genomes.

----------------
Trimming
----------------

Trimming is performed on both metaG and metaT reads. Trimming is performed by trimmomatic. Trimmomatic may produce
singletons from paired-end data, when one read is completely removed due to quality reasons. The singletons produced
from the first and second reads are concatenated into one file. After the trimmometic step, there are therefore always
three files for paired-end data, ``r1``, ``r2``, and ``se``.

A user-defined step is the removal of trailing Gs that are commonly introduced by Nextseq machines when the sequenced DNA
fragment is shorter than the number of bases added during the sequencing run. This is accomplished by cutadapt and is performed
on all reads, if requested in the config file.

------------
rRNA removal
------------

rRNA reads are separated from other metaT reads by sortmeRNA. The reason is that rRNA is highly abundant in total rRNA but doesn't
assemble readily with the default settings of the assemblers in IMP3. Very commonly, rRNA is actually depleted during library
preparation, with different success for different source organisms, making the rRNA abundance even less interpretable. While the rRNA
removal step removes the rRNA reads from futher processing, they are kept in a separate file for potential use outside of IMP3.

-------------------
Reference filtering
-------------------

Commonly, users will want to remove reads that map to one or more reference genomes, e.g. a host genome in a gut
microbiome or a known contaminant. IMP3 achieves this step for metaG and metaT reads by mapping against the files chosen
by the user with BWA. Only reads that do not map are kept. BWA is actually run independently on the paired-end data and singletons.
Both partners of a set of paired reads where one partner maps to the host genome are removed from the final data set.
If the user supplies more than one reference genome for filtering, the reads that did not map to the first reference will be
mapped against the second. The reads that did not map to this one will be mapped against the third reference and so on.
Currently reads that map to the reference genome(s) are not kept, nor are the alignments.

-----------
Cleaning up
-----------

Input and intermediary FASTQ files are gzipped at the end of the preprocessing step. 

