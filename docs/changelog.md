# Changelog

IMP3 is the successor of IMP (Narayanasamy et al., *Genome Biol* 2016;17:260). Code DOI for all versions:
[doi:10.5281/zenodo.24005282](https://doi.org/10.5281/zenodo.24005282); Apptainer image:
[doi:10.5281/zenodo.24003090](https://doi.org/10.5281/zenodo.24003090).

## unreleased (3.1.0)

- **Read filtering without false host removal** (`filtering.method`). `bwa_any`, the behaviour up to 3.0.3, removed a
  pair whenever BWA-MEM2 reported a mate as mapped; with short-insert libraries its pairing and mate rescue mark many
  microbial pairs as mapped: on 10 Expobiome stool metagenomes 12–31 % of the pairs, against 0.02–1.1 % human reads
  by Deacon and Hostile, with a GC-dependent loss (low-GC Bacillota most). New methods:
  - `deacon` (default of the configuration template): the host genome and transcriptome steps use
    [Deacon](https://github.com/bede/deacon) 0.20 and the human pangenome index `panhuman-1`; PhiX and other filter
    references use BWA-MEM2 with the rule of `bwa`;
  - `bwa`: BWA-MEM2 for all references, a read counts as aligned only when ≥ 80 % of it aligns (≥ 30 bases, ≤ 5 %
    mismatches), and the mate of a removed read is no longer lost (`scripts/host_filter_reads.py`).
  Both remove 0.01–0.1 % on these samples. Configurations without `filtering.method` keep `bwa_any`, so earlier runs
  can be reproduced. The installer fetches the Deacon index with `deacon: index: panhuman-1`; `deacon` is added to
  the IMP_mapping environment (the Apptainer image needs a rebuild).

- **Single-end metatranscriptomes** (single-omic mt runs with one read file): the preprocessing rules now use the
  same compressed file names as the paired-end rules (trimming, RiboDetector rRNA separation, host filtering and the
  final `mt.se.preprocessed.fq.gz` no longer disagree on `.fq` / `.fq.gz` and on step names), so a single-end stranded
  total-RNA library runs through preprocessing, assembly and annotation.
- **`featureCountsStranding: auto`**: the strand of a library is detected from its reads on the assembly (`-s 1` vs
  `-s 2` on the CDS; ≥ 3× decides, otherwise unstranded) and used by all featureCounts rules;
  `Stats/<omic>/<omic>.strandedness.json` keeps the counts.

- **Remapping step** (`steps: ... remap`, `remap: assembly_dir`): a run's preprocessed reads (e.g. an extracellular
  fraction) are mapped to the assembly of another IMP3 run of the same sample and counted on its annotation (genes,
  Mantis databases, per contig) and per MAG (reads, depth, Mantis terms), with strand detection (`stranding: auto`);
  outputs in `remapping/`. Replaces the manual Expobiome remapping scripts.

- **Term counts from the merged annotation**: the per-database counts (`<omic>.mantis.<db>.counts.tsv`, per contig,
  per MAG) now use the merged annotation `annotation_CDS_RNA_hmms.gff` (all sources), split into one line per gene and
  term; Pfam domains reported by Mantis under both accession and name (`PF00483` / `NTP_transferase`) are counted once,
  as accession (NAME / ACC of `Pfam-A.hmm`); KEGG pathways given as `ko` and `map` are counted once (`ko`). Before, a
  domain could be counted under two terms with different values, and with `-O --fraction` a read was split over the
  duplicates.

- **Bakta fork moved to `genomeanalysislcsb/bakta`** (from `ohickl/bakta`, branch `dynamic_db_fetch`), with a fix
  for assemblies without tmRNA: the fork read Aragorn's summary line by column and failed (`ValueError: 'found'`)
  when no chunk had a tmRNA, e.g. small single-end metaT assemblies. `annotate` also gives Bakta no more threads
  than the assembly has contigs.

- **`kraken_reads` for single-end metatranscriptomes**: the rRNA reads to add were a param, so they were missing
  from the shadow directory and the rule stopped without a message. They are now an input (`.fq.gz`).

## 3.0.3 (2026-11-23)

Fixes in the launcher and the summary statistics. The 3.0.0 Apptainer image is unchanged and still used.

- `imp3 run`: all bound folders are bound while probing the temporary root. Before, the probe of the second and
  later samples bound only the folders not yet checked and reported the temporary folder as not shared with the
  container, so a multi-sample `imp3 run` stopped after the first sample.
- Summary statistics (`summarize_all_stats.R`): the mapping rows are read from `contig_flagstat.txt` for both
  samtools layouts (16 lines since samtools 1.13, 13 lines before). With samtools ≥ 1.13 the fixed positions read the
  wrong lines (mapped segments 0, both / one partner mapped taken from the read1 / read2 counts). Rerun the summary
  rule to correct the statistics of existing runs.
- Docs: link to the documentation of the unpublished 2020 version (imp3.readthedocs.io); licence line names the
  LCSB; README and docs give the code DOI for all versions.

## 3.0.2 (2026-11-19)

- Summary statistics: bin quality tiers (HQ / MQ / LQ drafts) and the binned length were always 0, because the bins
  table now has the CheckM2 / MetaWRAP header (Completeness, Contamination, Genome_Size); columns are matched
  case-insensitively.

## 3.0.1 (2026-11-18)

Launcher fixes from the first large run (80 metagenomes on ULHPC). Uses the 3.0.0 image.

- `imp3 run --submit -j N`: at most N sample masters at a time, the others chained with `--dependency=afterany`
  (a failed sample does not hold up its chain); masters get `--time 24:00:00` in the default QoS.
- `imp3 init --site ulhpc --container`: `pre_commands` load the Apptainer module where `module` exists; container
  mount points are read from the image label, `tmp_root` is set to `$SCRATCH`, no site paths are built in.
- `mantis_annot` works in a job-private folder in the shadow directory instead of `/dev/shm`.
- Docs: Mantis out of memory (oom_kill after finishing) → run `mantis_annot` on large-memory nodes.

## 3.0.0 (2026-11-06)

First release of IMP3, with an Apptainer image on Zenodo.

- Preprocessing: fastp (adapter detection, quality trimming, pair merging); BWA-MEM2 against one or more
  contaminants and a host genome and transcriptome (default PhiX and human T2T-CHM13v2.0), batch size fixed so that
  results do not depend on the number of threads; RiboDetector for rRNA reads of the metatranscriptome; Nonpareil.
- Assembly: iterative metaG + metaT co-assembly with MEGAHIT (metaSPAdes and the long-read hybrid assembly removed).
- Annotation and quantification: Bakta; Mantis consensus annotation (KOfam, Pfam, eggNOG, TCDB and others);
  featureCounts per term and omic.
- Taxonomy: Kraken 2 with Bracken, mOTUs, contig taxonomy.
- Binning: CONCOCT, MetaBAT 2, MetaDecoder, SemiBin2, COMEBin and binny, combined by DAS Tool and refined with
  MetaWRAP; CheckM2 quality, GTDB-Tk taxonomy, MIMAG tiers, GRiD.
- `imp3` launcher: `init`, `run`, `status`, `tools`, `containerize`; one project with many samples on Slurm or a
  single machine, natively with conda or in an Apptainer image that holds all environments.
- Documentation renewed with MkDocs (installation, databases, inputs, workflow steps, outputs, FAQ, credits);
  references for every tool; gzipped test reads.
