
.. _step_assembly:

====================
IMP3 steps: Assembly
====================

IMP3 uses the iterative assembly approach of the original IMP. Like IMP, IMP3 can perform hybrid assemblies of metagenomic
and metatranscriptomic reads. In addition, IMP3 can perform purely metagenomic reads, ignoring metatranscriptomics reads. IMP3
can also perform another kind of hybrid assembly, namely of short and long metagenomic reads. Obviously, IMP3 also performs
iterative assemblies of only metagenomic or only metatranscriptomic reads.

After the assembly, reads are mapped back to the assembly.

.. _step_mono_iter_assembly:

-----------------------------------------------------------
Iterative assembly: metagenomic OR metatranscriptomic reads
-----------------------------------------------------------

In the simplest case, the user provides only one kind of short reads (metaG or metaT). IMP will then use the assembler defined
by the user (`Megahit <http://www.metagenomics.wiki/tools/assembly/megahit>`_  
or `MetaSpades <http://cab.spbu.ru/software/meta-spades/>`_) to try to assemble all reads. After the assembly, the reads are mapped
back to the assembled contigs. Reads that did not map will be given to the same assembler again.
The set of contigs from the second assembly will be merged with the first set using the overlap assembler
`CAP3 <http://seq.cs.iastate.edu/cap3.html>`_. The IMP developers
have `shown <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1116-8>`_ that this approach leads to longer
contigs without compromising assembly quality. Finally, the metaG or metaT are mapped
back to the final set of contigs by `BWA <http://bio-bwa.sourceforge.net/bwa.shtml>`_.

The same approach is performed for the metaG reads, if the user provides metaG and metaT but does not choose ``hybrid`` as assembly option.
In this case, both metaG and metaT reads are mapped to the final metagenomic assembly.

-------------------------------------------------
Hybrid assembly: long and short metagenomic reads
-------------------------------------------------

One of the assemblers built into IMP3, `MetaSpades <http://cab.spbu.ru/software/meta-spades/>`_,
is able to perform long/short-hybrid assemblies. If the user provides both kinds
of metagenomic reads, the two sets of reads are co-assembled in the first round of the iterative assembly.
Only the short reads are mapped back to the assembled contigs and, as described :ref:`above <step_mono_iter_assembly>`,
a second assembly will be attempted with the unmapped reads, and both sets of contigs will eventually be merged. 

.. _step_hybrid_iter_assembly:

------------------------------------------------------------------------------
Iterative multi-omic hybrid assembly: metagenomic and metatranscriptomic reads
------------------------------------------------------------------------------

This approach developed for the `original IMP <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1116-8>`_.
It is currently only implemented for the `Megahit <http://www.metagenomics.wiki/tools/assembly/megahit>`_ assembler.
First, the metatranscriptomic reads are assembled. The metatranscriptomic reads are mapped back to the contigs using
`BWA <http://bio-bwa.sourceforge.net/bwa.shtml>`_.
Reads that don't map are extracted and assembly is attempted a second time. The contigs from both assemblies are
supplied as reads to the hybrid assembly, together with all metagenomic and metatranscriptomic reads. All metagenomic and
metatranscriptomic reads are mapped against the resulting contigs. Metagenomic and metatranscriptomic reads that
don't map are extracted and are co-assembled. The resulting contigs are merged with the first set of hybrid contigs using
the `CAP3 <http://seq.cs.iastate.edu/cap3.html>`_ overlap assembler.

Finally, both metaG and metaT reads are mapped back to the final set of contigs.

