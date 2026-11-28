# IMP3: to do

Open points that are known but not scheduled. Add an item when something is removed, postponed or found
missing; remove it when it is done (the commit says so).

## Removed, could come back

- Removed together with their (never installed) environments: MaxBin and abawaca binners, BBNorm coverage
  normalization of metaT reads, EukDetect, HMMER annotation (`annotation: hmmer`, `hmm_DBs`), the MongoDB summary
  (`summary_steps: DB`), and GRiD (already disabled). Old rules are in the history before the commit that removed
  them; configs that ask for them stop with a message.

## Inputs

- metaG single-end reads: the Snakefile includes `modules/single_omics/mg/Preprocessing_SE.smk`, which does not
  exist (there is a half-finished `PreprocessingSE.smk` that asks for r1/r2 files); mg-only single-end runs stop.

## Assembly

- **metaSPAdes as assembler (removed in 04d0fa5).** Bring it back as an alternative to MEGAHIT, at least for
  metaG-only samples. The old rules are in the history (`git show 04d0fa5^:workflow/rules/Assembly/common/metaspades.smk`,
  `.../hybrid/metaspades.smk`). What did not work: the hybrid rules ran `METASPADES_ASSEMBLY_SHELL` as a command
  (the shell text was never filled in), and the integrative module includes `Assembly/hybrid/<assembler>.hybrid.smk`,
  which never existed for metaSPAdes. To do: rules for metaG (and, if useful, metaT with rnaSPAdes), the
  `assembler` value in the schema, cluster configs (big memory), docs; `spades` is still in `IMP_assembly.yaml`.
- **Long reads.** They were only used by the metaSPAdes hybrid assembly (`LongReads` now stops with a message).
  Options: metaSPAdes hybrid again, or Flye (already in `IMP_assembly.yaml`) followed by the short-read polishing
  and merging of the current workflow.

## Analysis

- **Fractional counting of overlapping terms.** All featureCounts rules use `-O --fraction -M --fraction`. For the
  term counts (`<omic>.mantis.<db>.counts.tsv`, per contig, per MAG, also in `remap`) this splits a read over all
  terms of one database on the same gene: a read on a gene with two Pfam domains gives 0.5 to each, a gene with one
  domain gives 1, so term counts depend on how many terms a gene carries. Option: `-O` without `--fraction` for the
  term counts (each term gets the whole read; multi-mapping stays fractional), keep `-O --fraction` for the gene
  counts. Left as is in 3.1.0 (the duplicate Pfam name / accession terms that made it visible are merged since
  then); decide before the next comparison of term counts across runs.
- **imp-explorer and Pfam accessions.** Since 3.1.0 the Pfam term counts use accessions only (`PF00483`, no names);
  check that imp-explorer's Pfam layers and term names read them (and old releases with names).

## Installation and tests

- `config/config_test.yaml` is older than the workflow: it lacks keys the rules read (e.g. `mantis: params`), so the
  test run stops while parsing. Regenerate it from `config/config.imp.yaml` with the test reads and test paths.
- Bakta copies its database (72 GB on ULHPC) to `/dev/shm` (`--tmp_db_path`), which counts against the job's memory;
  with `max_threads: 16` (test runs) the annotate job gets too little memory and is OOM-killed while Bakta seems to
  hang ("predict CRISPR arrays..." is the last log line). Size the annotate job's memory from the database size (or
  read the database from disk when the job is small), and say so in the docs of `max_threads`.
- `localrules` names `prepare_input_data`, which no longer exists (Snakemake warns at every start).
- Snakemake 7 (`imp_sm_classic`) with mamba >= 2: Snakemake creates the environment folder before
  `mamba env create --prefix`, and mamba 2 aborts ("Non-conda folder exists at prefix"). Workaround: create the
  environment with `micromamba create -p <conda-prefix>/<hash>_ -f <conda-prefix>/<hash>_.yaml`, or use
  `--conda-frontend conda` / mamba < 2 in the Snakemake environment. Fix in the installation docs (or move to
  Snakemake 8).

## Container

- Tool caches in the sample folder: in the container Snakemake sets `--home` to the sample's output folder, so
  matplotlib, Java, Keras and GNU parallel leave `.cache/`, `.java/`, `.keras/`, `.parallel/` there (harmless, a few
  files). Could point HOME to the run's temporary folder instead.
  The container run also left `nltk_data/` (NLTK stopwords downloaded by a tool) in the sample folder.

## smORFs

- Metaproteomics search databases with the smORFs (the paper's MP_DB: Bakta proteins + smORFs + human RefSeq + cRAP,
  deduplicated, one ID per identical protein across samples) were part of the former Snakemake version of the
  smORF module (git history: workflow/rules/smORF/proteomics_db.smk) and are not in `imp3 smorf` yet.
- `flank` > 0 (~60 nt) to recover smORFs in gaps that touch a neighbouring feature (690 on BHC01 found only by the
  paper's whole-assembly gene calling); compare on a few samples and choose the default.
- An IMP3 step `smorf` (rule calling `imp3_smorf.py` before the final compression) so new runs can include it.
- `IMP_smorf` environment into the next Apptainer image.
