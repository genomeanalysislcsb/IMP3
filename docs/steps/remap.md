# Remapping to another run

A study can sequence further fractions of the same samples, e.g. the extracellular DNA and RNA next to the
intracellular material of the main runs. With `remap` in `steps`, a run maps its preprocessed reads to the assembly
of another IMP3 run of the same sample and counts them on that run's annotation and MAGs, with the same featureCounts
settings as the main workflow. The fraction is then directly comparable with the other run's own reads, gene by gene,
term by term and MAG by MAG.

```yaml
steps: "preprocessing remap"            # e.g. an extracellular metaT library, reads only
remap:
  assembly_dir: /path/IMP/intra_MGMT/S01
  stranding:
    mt: auto
```

A fraction can also get a full run of its own (`steps: "preprocessing assembly analysis taxonomy binning remap"`):
its own assembly and MAGs, plus the remapping.

| step | output |
|---|---|
| bwa-mem2 index of the other run's assembly (`remap.assembly`) | `remapping/Assembly/remap.assembly.fa.*` |
| mapping (pairs and singletons, or single-end reads) | `remapping/Assembly/<omic>.reads.sorted.bam`, `remapping/Stats/<omic>/<omic>.flagstat.txt` |
| reads and depth per contig (`samtools coverage`) | `remapping/Stats/<omic>/<omic>.contig_coverage.tsv` |
| strand of the library (`stranding: auto`: `-s 1` vs `-s 2` on the CDS, ≥ 3× decides, else 0) | `remapping/Stats/<omic>/<omic>.strandedness.json` |
| genes (`ID`, `Name` × CDS, rRNA, tRNA, ncRNA, CRISPR, …) | `remapping/annotation/<omic>.<ID\|Name>.<feature>_counts.tsv` |
| Mantis terms per database | `remapping/annotation/<omic>.mantis.<db>.counts.tsv` |
| the same per contig | `remapping/annotation/per_contig_counts/` |
| reads per MAG (`remap.bins`) | `remapping/Stats/<omic>/<omic>.bin_abundance.tsv` (contigs, length, reads, share, mean depth; `unbinned`) |
| Mantis terms per MAG | `remapping/annotation/per_bin_counts/<omic>.mantis.<db>.per_bin.counts.tsv` (bin, term, count) |

The Mantis counts use the other run's per-database GFF (`Analysis/annotation/intermediary/mantis/ref_gff/`) when it
is still there, as its own counts did, and otherwise the merged annotation `annotation_CDS_RNA_hmms.gff` with the
database as attribute. featureCounts counts multi-mapping and overlapping reads fractionally (`-O --fraction -M
--fraction`), paired reads with `-p`, as in the main workflow. imp-explorer reads the counts with
`imp-explorer fraction --imp-root <folder of the fraction runs> --suffix ext` (omics `<omic>_ext_remap`).
