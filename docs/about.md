# History and credits

- **IMP** (2016) introduced the reproducible, integrated analysis of coupled metaG and metaT data with an
  iterative multi-omic hybrid assembly:
  Narayanasamy S, Jarosz Y, Muller EEL, Heintz-Buschart A, *et al.* IMP: a pipeline for reproducible
  reference-independent integrated metagenomic and metatranscriptomic analyses.
  *Genome Biology* 17, 260 (2016). <https://doi.org/10.1186/s13059-016-1116-8>
- **IMP3** re-implemented IMP in Snakemake with conda environments and extended it with binning, read- and
  contig-based taxonomy, long reads and modular steps (first documentation 2020, Anna Heintz-Buschart).
- **IMP3 3.0** (current) replaced Trimmomatic by fastp and Prokka by Bakta, added Mantis annotation, several
  binners with DAS Tool / MetaWrap refinement, CheckM2, GTDB-Tk, GRiD, mOTUs, Nonpareil and per-sample
  statistics tables.

The development history moved from the LCSB GitLab to
[GitHub](https://github.com/genomeanalysislcsb/IMP3) with all commits and authors.

## Licence

MIT licence, © Helmholtz Centre for Environmental Research (UFZ) & Luxembourg Centre for Systems Biomedicine
(LCSB).

## Citing the tools

When you publish IMP3 results, cite IMP and the tools of the steps you ran. The main ones:
fastp, BWA-MEM2, SortMeRNA, MEGAHIT, quickmerge / CAP3, Bakta, Mantis, HMMER, featureCounts,
Kraken2, Bracken, mOTUs, GTDB-Tk, CheckM2, DAS Tool, MetaWrap and the binners you used (CONCOCT, MetaBAT 2,
MetaDecoder, SemiBin, COMEBin, binny, MaxBin).
