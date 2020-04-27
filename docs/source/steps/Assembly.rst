
.. _step_assembly:

====================
IMP3 steps: Assembly
====================

IMP3 uses the iterative assembly approach of the original IMP. Like IMP, IMP3 can perform hybrid assemblies of metagenomic
and metatranscriptomic reads. In addition, IMP3 can perform purely metagenomic reads, ignoring metatranscriptomics reads. IMP3
can also perform another kind of hybrid assembly, namely of short and long metagenomic reads. Obviously, IMP3 also performs
iterative assemblies of only metagenomic or only metatranscriptomic reads.

After the assembly, reads are mapped back to the assembly.

-----------------------------------------------------------
Iterative assembly: metagenomic OR metatranscriptomic reads
-----------------------------------------------------------

In the simplest case, the user provides only one kind of short reads (metaG or metaT). IMP will then use the assembler defined
by the user to assemble all reads. After the assembly, the reads are mapped back to the assembled contigs. Reads that did not map
will be assembled. ??? The second set of contigs will be merged with the first set using the overlap assembler CAP3. The IMP developers
have shown that this approach leads to longer contigs without compromising assembly quality. Finally, the metaG or metaT are mapped
back to the final set of contigs by BWA.

The same approach is performed for the metaG reads, if the user provides metaG and metaT but chooses nonhybrid as assembly option.
In this case, metaG and metaT reads are mapped to the final metagenomic assembly by BWA.

-------------------------------------------------
Hybrid assembly: long and short metagenomic reads
-------------------------------------------------

One of the assemblers built into IMP3, metaspades, is able to perform long/short-hybrid assemblies. If the user provides both kinds
of metagenomic reads, the two sets of reads are co-assembled in the first round of the iterative assembly ???
Only the short reads are mapped back to the assembled contigs ???

------------------------------------------------------------------------------
Iterative multi-omic hybrid assembly: metagenomic and metatranscriptomic reads
------------------------------------------------------------------------------

This approach developed for the original IMP. First, the metatranscriptomic reads are assembled using megahit. The resulting reads...

Finally, both metaG and metaT reads are mapped back to the final set of contigs.

