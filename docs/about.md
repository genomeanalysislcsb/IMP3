# History and credits

- **IMP** (2016) introduced the reproducible, integrated analysis of coupled metaG and metaT data with an
  iterative multi-omic hybrid assembly:
  Narayanasamy S, Jarosz Y, Muller EEL, Heintz-Buschart A, *et al.* IMP: a pipeline for reproducible
  reference-independent integrated metagenomic and metatranscriptomic analyses.
  *Genome Biology* 17, 260 (2016). <https://doi.org/10.1186/s13059-016-1116-8>
- **IMP3** re-implemented IMP in Snakemake with conda environments and extended it with binning, read- and
  contig-based taxonomy, long reads and modular steps (first documentation 2020, Anna Heintz-Buschart).
- **IMP3 3.0** (current) replaced Trimmomatic by fastp and Prokka by Bakta, added Mantis annotation, several
  binners with DAS Tool / MetaWrap refinement, CheckM2, GTDB-Tk, mOTUs, Nonpareil and per-sample
  statistics tables.

The development history moved from the LCSB GitLab to
[GitHub](https://github.com/genomeanalysislcsb/IMP3) with all commits and authors.

## Licence

MIT licence, © Helmholtz Centre for Environmental Research (UFZ) & Luxembourg Centre for Systems Biomedicine
(LCSB).

## Citing the tools

When you publish IMP3 results, cite IMP, Snakemake and the tools of the steps you ran (the binners only if they
were in `binning.binners`). The versions are in the environment files (`workflow/envs/*.yaml`).

### Workflow

| tool | reference |
|---|---|
| IMP | Narayanasamy S *et al.* IMP: a pipeline for reproducible reference-independent integrated metagenomic and metatranscriptomic analyses. *Genome Biology* 17, 260 (2016). [doi:10.1186/s13059-016-1116-8](https://doi.org/10.1186/s13059-016-1116-8) |
| Snakemake | Mölder F *et al.* Sustainable data analysis with Snakemake. *F1000Research* 10, 33 (2021). [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2) |
| Apptainer / Singularity (container runs) | Kurtzer GM, Sochat V, Bauer MW. Singularity: scientific containers for mobility of compute. *PLoS ONE* 12, e0177459 (2017). [doi:10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459) |

### Preprocessing

| tool | reference |
|---|---|
| fastp | Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34, i884–i890 (2018). [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560) |
| BWA-MEM2 | Vasimuddin M, Misra S, Li H, Aluru S. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. *IEEE IPDPS 2019*, 314–324. [doi:10.1109/IPDPS.2019.00041](https://doi.org/10.1109/IPDPS.2019.00041) |
| SAMtools, BCFtools | Danecek P *et al.* Twelve years of SAMtools and BCFtools. *GigaScience* 10, giab008 (2021). [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| T2T-CHM13 (human reference for host filtering) | Nurk S *et al.* The complete sequence of a human genome. *Science* 376, 44–53 (2022). [doi:10.1126/science.abj6987](https://doi.org/10.1126/science.abj6987) |
| RiboDetector (rRNA reads, metaT) | Deng Z-L, Münch PC, Mreches R, McHardy AC. Rapid and accurate identification of ribosomal RNA sequences via deep learning. *Nucleic Acids Research* 50, e60 (2022). [doi:10.1093/nar/gkac112](https://doi.org/10.1093/nar/gkac112) |
| FastQC | Andrews S. FastQC: a quality control tool for high throughput sequence data (2010). <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/> |
| Nonpareil | Rodriguez-R LM, Gunturu S, Tiedje JM, Cole JR, Konstantinidis KT. Nonpareil 3: fast estimation of metagenomic coverage and sequence diversity. *mSystems* 3, e00039-18 (2018). [doi:10.1128/mSystems.00039-18](https://doi.org/10.1128/mSystems.00039-18) |
| SeqKit | Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLoS ONE* 11, e0163962 (2016). [doi:10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962) |

### Assembly

| tool | reference |
|---|---|
| MEGAHIT | Li D, Liu C-M, Luo R, Sadakane K, Lam T-W. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* 31, 1674–1676 (2015). [doi:10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033) |
| quickmerge (`assembly.merge: quickmerge`) | Chakraborty M, Baldwin-Brown JG, Long AD, Emerson JJ. Contiguous and accurate *de novo* assembly of metazoan genomes with modest long read coverage. *Nucleic Acids Research* 44, e147 (2016). [doi:10.1093/nar/gkw654](https://doi.org/10.1093/nar/gkw654) |
| MUMmer4 (used by quickmerge) | Marçais G *et al.* MUMmer4: a fast and versatile genome alignment system. *PLoS Computational Biology* 14, e1005944 (2018). [doi:10.1371/journal.pcbi.1005944](https://doi.org/10.1371/journal.pcbi.1005944) |
| CAP3 (`assembly.merge: cap3`) | Huang X, Madan A. CAP3: a DNA sequence assembly program. *Genome Research* 9, 868–877 (1999). [doi:10.1101/gr.9.9.868](https://doi.org/10.1101/gr.9.9.868) |
| Flye (`assembly.merge: flye`) | Kolmogorov M, Yuan J, Lin Y, Pevzner PA. Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology* 37, 540–546 (2019). [doi:10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8) |
| CD-HIT | Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics* 28, 3150–3152 (2012). [doi:10.1093/bioinformatics/bts565](https://doi.org/10.1093/bioinformatics/bts565) |
| mosdepth | Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. *Bioinformatics* 34, 867–868 (2018). [doi:10.1093/bioinformatics/btx699](https://doi.org/10.1093/bioinformatics/btx699) |
| hts-nim (hts-nim-tools) | Pedersen BS, Quinlan AR. hts-nim: scripting high-performance genomic analyses. *Bioinformatics* 34, 3387–3389 (2018). [doi:10.1093/bioinformatics/bty358](https://doi.org/10.1093/bioinformatics/bty358) |
| BEDTools | Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26, 841–842 (2010). [doi:10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033) |

### Annotation and counts

| tool | reference |
|---|---|
| Bakta | Schwengers O *et al.* Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics* 7, 000685 (2021). [doi:10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685) |
| Mantis | Queirós P, Delogu F, Hickl O, May P, Wilmes P. Mantis: flexible and consensus-driven genome annotation. *GigaScience* 10, giab042 (2021). [doi:10.1093/gigascience/giab042](https://doi.org/10.1093/gigascience/giab042) |
| KOfam (Mantis reference) | Aramaki T *et al.* KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. *Bioinformatics* 36, 2251–2252 (2020). [doi:10.1093/bioinformatics/btz859](https://doi.org/10.1093/bioinformatics/btz859) |
| Pfam (Mantis reference) | Mistry J *et al.* Pfam: the protein families database in 2021. *Nucleic Acids Research* 49, D412–D419 (2021). [doi:10.1093/nar/gkaa913](https://doi.org/10.1093/nar/gkaa913) |
| eggNOG (Mantis reference) | Huerta-Cepas J *et al.* eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource. *Nucleic Acids Research* 47, D309–D314 (2019). [doi:10.1093/nar/gky1085](https://doi.org/10.1093/nar/gky1085) |
| TCDB (Mantis reference) | Saier MH Jr *et al.* The Transporter Classification Database (TCDB): 2021 update. *Nucleic Acids Research* 49, D461–D467 (2021). [doi:10.1093/nar/gkaa1004](https://doi.org/10.1093/nar/gkaa1004) |
| featureCounts | Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30, 923–930 (2014). [doi:10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656) |

### Taxonomy

| tool | reference |
|---|---|
| Kraken2 | Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biology* 20, 257 (2019). [doi:10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0) |
| Bracken | Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science* 3, e104 (2017). [doi:10.7717/peerj-cs.104](https://doi.org/10.7717/peerj-cs.104) |
| KrakenTools | Lu J *et al.* Metagenome analysis using the Kraken software suite. *Nature Protocols* 17, 2815–2839 (2022). [doi:10.1038/s41596-022-00738-y](https://doi.org/10.1038/s41596-022-00738-y) |
| mOTUs | Ruscheweyh H-J *et al.* Cultivation-independent genomes greatly expand taxonomic-profiling capabilities of mOTUs across various environments. *Microbiome* 10, 212 (2022). [doi:10.1186/s40168-022-01410-z](https://doi.org/10.1186/s40168-022-01410-z) |
| GTDB-Tk | Chaumeil P-A, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk v2: memory friendly classification with the genome taxonomy database. *Bioinformatics* 38, 5315–5316 (2022). [doi:10.1093/bioinformatics/btac672](https://doi.org/10.1093/bioinformatics/btac672) |
| GTDB | Parks DH *et al.* GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research* 50, D785–D794 (2022). [doi:10.1093/nar/gkab776](https://doi.org/10.1093/nar/gkab776) |

### Binning

| tool | reference |
|---|---|
| CONCOCT | Alneberg J *et al.* Binning metagenomic contigs by coverage and composition. *Nature Methods* 11, 1144–1146 (2014). [doi:10.1038/nmeth.3103](https://doi.org/10.1038/nmeth.3103) |
| MetaBAT 2 | Kang DD *et al.* MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ* 7, e7359 (2019). [doi:10.7717/peerj.7359](https://doi.org/10.7717/peerj.7359) |
| MetaDecoder | Liu C-C *et al.* MetaDecoder: a novel method for clustering metagenomic contigs. *Microbiome* 10, 46 (2022). [doi:10.1186/s40168-022-01237-8](https://doi.org/10.1186/s40168-022-01237-8) |
| SemiBin2 | Pan S, Zhao X-M, Coelho LP. SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. *Bioinformatics* 39, i21–i29 (2023). [doi:10.1093/bioinformatics/btad209](https://doi.org/10.1093/bioinformatics/btad209) |
| COMEBin | Wang Z, You R, Han H, Liu W, Sun F, Zhu S. Effective binning of metagenomic contigs using contrastive multi-view representation learning. *Nature Communications* 15, 585 (2024). [doi:10.1038/s41467-023-44290-z](https://doi.org/10.1038/s41467-023-44290-z) |
| binny | Hickl O, Queirós P, Wilmes P, May P, Heintz-Buschart A. binny: an automated binning algorithm to recover high-quality genomes from complex metagenomic datasets. *Briefings in Bioinformatics* 23, bbac431 (2022). [doi:10.1093/bib/bbac431](https://doi.org/10.1093/bib/bbac431) |
| DAS Tool | Sieber CMK *et al.* Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nature Microbiology* 3, 836–843 (2018). [doi:10.1038/s41564-018-0171-1](https://doi.org/10.1038/s41564-018-0171-1) |
| MetaWRAP | Uritskiy GV, DiRuggiero J, Taylor J. MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis. *Microbiome* 6, 158 (2018). [doi:10.1186/s40168-018-0541-1](https://doi.org/10.1186/s40168-018-0541-1) |
| CheckM2 | Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods* 20, 1203–1212 (2023). [doi:10.1038/s41592-023-01940-w](https://doi.org/10.1038/s41592-023-01940-w) |
