
.. _output_analysis:

=====================
IMP3 output: Analysis
=====================

----------
Annotation
----------

The most important results of the **analysis** step are in the **annotation** subdirectory.

``annotation_CDS_RNA_hmms.gff``
contains all gene annotations (including mRNAs, rRNAs and tRNAs) with the functional annotations. The file is in `standard
.gff format <https://www.ensembl.org/info/website/upload/gff.html>`_, this means nine columns (1: contig ID, 2: source
software, 3: kind of feature ("CDS" for protein-coding regions,
"rRNA" for rRNAs, "tRNA" for tRNAs, and "tmRNA" for tRNAs with coding regions), 4: left-most coordinate of feature (1-based), 5:
right-most coordinate of the feature (1-based, inclusive), 6: a score ("." if no score is reported), 7: the
feature's direction or sense ("+" or "-"), 8: frame (0 for CDS, "." for other features), 9:
attributes (each attribute starts with a key followed by "=", e.g. "ID=", attributes are separated by ";" - the
most interesting attributes are **ID**, which is used in other gene-based outputs, **partial**, which contains information from the
gene-caller, prodigal, about whether the open reading frames are complete (00 means both start and stop codons were found, 11 means
neither were found, 01 means the right-most end is incomplete - i.e. missing stop-codon for +strand features, missing start-codon for -strand
features, 10 means the left-most end is incomplete - i.e. missing start-codon for +strand features, missing stop-codon for -strand
features), and **the results of the HMM searches**, named with the database name you provided in the config file, e.g. **essential=** for
the essential genes)).

The annotated features are used to determine the numbers of reads mapping to the features and to groups of features that share the
same functional annotation. The results are in same functional annotation. The results are in ``<mg|mt>.<CDS|rRNA>_counts.<tsv>``
for the single features and ``<mg|mt>.<database>.counts.<tsv>`` for the functional groups. The format of these files is unchanged
`featureCounts <http://bioinf.wehi.edu.au/featureCounts/>`_ output, i.e. there's a line starting with "#" that contains the
original call, followed by seven tab-separated columns with a header (1: feature or group ID, 2: contig or
contigs, separated by ";", 3: left-most coordinate (as described for the .gff format) or coordinates (separated by ";"), 4:
right-most coordinate (as described for the .gff format) or coordinates (separated by ";"), 5: feature direction or directions
(separated by ";"), 6: total length of feature(s), 7: total counted reads).
In addition, summaries of the numbers of
reads that were mapped and overlapped with the features are found in the respective files ending on ``tsv.summary``. 

The FASTA file for proteomics is called ``proteomics.final.faa``. The header for the protein predictions is in a format that
should work with most proteomics search engines and analysis tools, in particular the `MetaProteomeAnalyzer <http://www.mpa.ovgu.de/>`_ .
If you've chosen to add host proteins, you need to make sure your host protein file also conforms to the requirements.
An intermediary step (without host proteins) is in ``proteomics.proteins.faa``.

The proteomics file is a cleaned version of ``prokka.faa``. In addition to the amino-acid
sequences, `prokka <https://github.com/tseemann/prokka>`_ also outputs
the CDS as ``prokka.ffn``, the contigs in ``prokka.fna`` and - with a slightly different header - in ``prokka.fsa``, as well
as a logfile ``prokka.log``, a summary of the number of analysed contigs and annotated features ``prokka.txt``, a tabular output
with all features ``prokka.tsv``, and a sort-of flat version of the same information in ``prokka.tbl``, and a gff file with
lots of commented lines (starting with "#") ``prokka.gff``, which is actually the foundation of the gff file described above.
An intermediary step between this file and the final gff is ``annotation.filt.gff`` which contains all the information of the
original prokka output minus the comments. Depending on your further steps, you may also see indices for some of the sequences
(``prokka.<faa|ffn>.<suffix>``). If you ran **binning**, the binning step will add a file containing the link between the gene IDs
as given by prokka and the format required by DASTool, called ``annotation_CDS_RNA_hmms.contig2ID.tsv``.

The functional annotation of the genes by `HMMer <http://hmmer.org/>`_ actually produces a number of intermediary outputs. With
the results summarized in ``annotation_CDS_RNA_hmms.gff``, the intermediary files are archived and compressed in
``intermediary.tar.gz``.

----
SNPs
----

IMP3 uses `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ to record positions where mapping reads
differ from the assembly consensus. The information is stored in the **snps** subfolder of the **Analysis** directory,
in ``<mg|mt>.variants.samtools.vcf.gz`` with the index file ``<mg|mt>.variants.samtools.vcf.gz.tbi``.


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

In addition to the annotation and SNP calling, the analysis step will also collect some stats for you and save them to the ``Stats
directory. 

