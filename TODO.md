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

## Container

- Tool caches in the sample folder: in the container Snakemake sets `--home` to the sample's output folder, so
  matplotlib, Java, Keras and GNU parallel leave `.cache/`, `.java/`, `.keras/`, `.parallel/` there (harmless, a few
  files). Could point HOME to the run's temporary folder instead.
  The container run also left `nltk_data/` (NLTK stopwords downloaded by a tool) in the sample folder.
