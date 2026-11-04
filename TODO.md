# IMP3: to do

Open points that are known but not scheduled. Add an item when something is removed, postponed or found
missing; remove it when it is done (the commit says so).

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

- Rules that `imp3 tools` reports as missing in the image but that a normal run does not need: `hmmsearch`
  (`annotation: hmmer`, not in `IMP_annotation.yaml`), `java` (VizBin), and the database setup rules
  (`bwa-mem2`, `prokka`, `mongod`, `checkm2`). Either add them to their environments / the image or remove the
  rules that are no longer used.
