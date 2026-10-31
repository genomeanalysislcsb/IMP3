# Inputs

IMP3 usually starts from raw paired-end reads, but it can also start from preprocessed reads, from an existing
assembly with alignments, or from an annotated assembly. All inputs are set under `raws` in the
[config file](configuration.md); paths should be absolute.

## Raw reads (default)

Two FASTQ files (gzipped or not) per omic, forward and reverse, reads in the same order:

```yaml
raws:
  Metagenomics: "/path/mg.r1.fastq.gz /path/mg.r2.fastq.gz"
  Metatranscriptomics: "/path/mt.r1.fastq.gz /path/mt.r2.fastq.gz"
  LongReads: ""
  LongReadTech: ""
  Contigs: ""
  Alignment_metagenomics: ""
  Alignment_metatranscriptomics: ""
  Gff: ""
```

Leave `Metatranscriptomics` (or `Metagenomics`) empty for a single-omic run. Single-end data: give one file per omic.

## Preprocessed reads

Give **three** files per omic, forward, reverse and singletons (last), and leave `preprocessing` out of `steps`:

```yaml
raws:
  Metagenomics: "/path/mg.r1.fq.gz /path/mg.r2.fq.gz /path/mg.se.fq.gz"
  Metatranscriptomics: "/path/mt.r1.fq.gz /path/mt.r2.fq.gz /path/mt.se.fq.gz"
```

## Long reads

Preprocessed long metaG reads, used together with short metaG reads in a metaSPAdes hybrid assembly
(`assembly.assembler: metaspades`):

```yaml
raws:
  Metagenomics: "/path/mg.r1.fastq.gz /path/mg.r2.fastq.gz"
  LongReads: "/path/longreads.fastq"
  LongReadTech: "nanopore"        # or pacbio
```

## Existing assembly

To annotate and bin an assembly made elsewhere, give the contigs (FASTA headers = contig names only, no spaces)
and either reads or alignments (BAM sorted by coordinate, same contig names), and leave `preprocessing` and
`assembly` out of `steps`:

```yaml
raws:
  Contigs: "/path/assembly.fasta"
  Alignment_metagenomics: "/path/mg.sorted.bam"
  Alignment_metatranscriptomics: "/path/mt.sorted.bam"
```

Add `Gff` (an annotation of these contigs) to run only binning.

## Steps and inputs

| steps | needs |
|---|---|
| all (`preprocessing assembly analysis binning taxonomy summary`) | raw reads |
| `assembly …` | preprocessed reads (three files per omic) |
| `analysis binning taxonomy …` | `Contigs` + reads or alignments |
| `binning …` | `Contigs`, alignments, `Gff` |
