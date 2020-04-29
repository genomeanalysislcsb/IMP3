.. _output_analysis:

=====================
IMP3 output: Analysis
=====================

----------
Annotation
----------

IMP3 writes the results of the :ref:`Analysis <step_analysis>` into the ``Analysis`` directory within in the defined ``outputdir``, see
:ref:`configuration <configuration>`. 

The gene annotation results are in the ``annotation`` subdirectory.

............
The GFF file
............

The `GFF3 <http://gmod.org/wiki/GFF3>`_ file ``annotation_CDS_RNA_hmms.gff`` is the final annotation file and contains all gene annotations (including mRNAs, rRNAs, tRNAs and CRISPR arrays) with the functional annotations and gene descriptions. 
The file is in `standard GFF format <https://www.ensembl.org/info/website/upload/gff.html>`_ with nine columns:

- 1: contig ID, 
- 2: source software, 
- 3: kind of feature ("CDS" for protein-coding regions, "rRNA" for rRNAs, "tRNA" for tRNAs, and "tmRNA" for tRNAs with coding regions), 
- 4: left-most coordinate of feature (1-based), 
- 5: right-most coordinate of the feature (1-based, inclusive), 
- 6: a score ("." if no score is reported), 
- 7: the feature's direction or sense ("+" or "-"), 
- 8: frame (0 for CDS, "." for other features), 
- 9: attributes (each attribute starts with a key followed by "=", e.g. "ID=", attributes are separated by ";" )
    
The most important attributes of column *9* are:

- **ID**, which is used in other gene-based outputs,
- **partial**, which contains information from the gene-caller (`Prodigal <https://github.com/hyattpd/Prodigal>`_) about whether the open reading frames are complete:
    - *00* means both start and stop codons were found, 
    - *11* means neither were found, 
    - *01* means the right-most end is incomplete - i.e. missing stop-codon for +strand features, missing start-codon for -strand features, 
    - *10* means the left-most end is incomplete - i.e. missing start-codon for +strand features, missing stop-codon for -strand features, and 
- **the results of the HMM searches**, named with the HMM database name provided in the **config file**, e.g. *essential=* for the essential genes.

The functional annotation of the genes by `HMMer <http://hmmer.org/>`_ produces a number of intermediary outputs. With the results summarized 
in ``annotation_CDS_RNA_hmms.gff``, the intermediary files are archived and compressed in ``intermediary.tar.gz``.
......................
Reads per gene / group
......................

The annotated features are used to determine the numbers of reads mapping to the features and to groups of features that share the
same functional annotation using `featureCounts <http://bioinf.wehi.edu.au/featureCounts/>`_:

- The results are written to ``<mg|mt>.<CDS|rRNA>_counts.<tsv>`` for single features and to ``<mg|mt>.<database>.counts.<tsv>`` for the functional groups. 
- The format of these files is unchanged `featureCounts <http://bioinf.wehi.edu.au/featureCounts/>`_ output, i.e. there's a line starting with "#" that contains the original call, followed by seven tab-separated columns with a header: 
    - 1: feature or group ID, 
    - 2: contig or contigs, separated by ";", 
    - 3: left-most coordinate (as described for the .gff format) or coordinates (separated by ";"), 
    - 4: right-most coordinate (as described for the .gff format) or coordinates (separated by ";"), 
    - 5: feature direction or directions (separated by ";"), 
    - 6: total length of feature(s), 7: total counted reads).
- In addition, summaries of the numbers of reads that were mapped and overlapped with the features are found in the respective files ending on ``tsv.summary``. 
- If only **metaT** reads were used as input, there will be no data for rRNA, because rRNA has been filtered out.

.......................................
Proteomics databases and gene sequences
.......................................

IMP3 output can be used for downstream **metaP** analysis:

- ``proteomics.final.faa``: the FASTA file for downstream **metaP** analysis. The FASTA header for the protein sequences is in a format that
should work with most proteomics search engines and analysis tools, in particular the `MetaProteomeAnalyzer <http://www.mpa.ovgu.de/>`_ .
- As an intermediary step (without host proteins), a file named ``proteomics.proteins.faa`` is generated.

The proteomics file is a cleaned version of ``prokka.faa``. `prokka <https://github.com/tseemann/prokka>`_ outputs a few more files:
- ``prokka.ffn`` with the CDS,
- ``prokka.fna`` with the contigs (also present in ``prokka.fsa`` with a slightly different header),
- ``prokka.log``, a logfile, 
- ``prokka.txt``, a summary of the number of analysed contigs and annotated features,
- ``prokka.tsv``, a tabular output with all features,
- ``prokka.tbl``, a sort-of flat version of the same information,
- ``prokka.gff``, a `GFF <http://gmod.org/wiki/GFF3>`_ file with lots of commented lines (starting with "#"), which is actually the foundation of the `GFF <http://gmod.org/wiki/GFF3>`_ file described above.

An intermediary step between this file and the final `GFF <http://gmod.org/wiki/GFF3>`_ is ``annotation.filt.gff`` which contains all the information of the
original prokka output minus the comments. Depending on the planned further analysis steps, the user may also see indices for some of the sequences
(``prokka.<faa|ffn>.<suffix>``). If the IMP3 :ref:`Binning <step_binning>` is run, IMP3 will add a file containing the link between the gene IDs
as given by `prokka <https://github.com/tseemann/prokka>`_ and the format required by `DASTool <https://github.com/cmks/DAS_Tool>`_, called ``annotation_CDS_RNA_hmms.contig2ID.tsv``.

----
SNPs
----

IMP3 uses `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ to call variants. `Samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ 
records positions where mapping reads differ from the assembly consensus. The information is stored in the ``snps`` subfolder of the ``Analysis`` 
directory in `VCF <https://en.wikipedia.org/wiki/Variant_Call_Format>`_ and `tabix <http://www.htslib.org/doc/tabix.html>`_ formats: ``<mg|mt>.variants.samtools.vcf.gz`` with the index file ``<mg|mt>.variants.samtools.vcf.gz.tbi``.


--------
Taxonomy
--------

The ``Analysis`` directory contains the output of the :ref:`taxonomy <output_taxonomy>` step. The exception is the output of `GTDBtk <https://github.com/Ecogenomics/GTDBTk>`_, 
which is run on the bins selected by `DASTool <https://github.com/cmks/DAS_Tool>`_. See the :ref:`next <output_taxonomy>` section for details.

.. _output_analysis_stats:

----------------------------
Stats from the Analysis step
----------------------------

In addition to the annotation and SNP calling, the analysis step will also collect some stats for you and save them to the **Stats** directory. 

