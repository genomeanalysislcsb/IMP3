# Analysis

Genes and RNAs are called on the contigs, annotated with functions, and the reads of each omic are counted per
gene and per functional term. Output: `Analysis/annotation/`, `Analysis/snps/`, `Stats/{mg,mt}/`.

## Gene calling (Bakta)

[Bakta](https://github.com/oschwengers/bakta) (an IMP3 fork, set up on first run) calls protein-coding genes,
including partial genes at contig ends, rRNA, tRNA, tmRNA, ncRNA and regulatory regions, CRISPR arrays, and oriC /
oriT sites. Bakta's own functional annotation is kept in `Analysis/annotation/bakta/`; IMP3 uses the filtered gene
set (`annotation.filt.gff`, `annotation.faa`) for the functional annotation below.

## Functional annotation (Mantis)

[Mantis](https://github.com/PedroQueirosLab/mantis) annotates every protein against all configured reference
databases and builds a consensus: KEGG orthologs, modules, pathways, BRITE, reactions and reaction classes, EC
numbers, Pfam, TIGRFAM, COG / arCOG, eggNOG, GO, CAZy, TCDB and essential genes (`mantis` in the
[config](../configuration.md#annotation-and-functions)). All annotations are merged into `Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz`.

## Read counts

[featureCounts](https://subread.sourceforge.net/) counts the reads of each omic (metaG unstranded, metaT
reverse-stranded: `featureCountsStranding`) on the merged annotation, with multi-mapping and overlapping reads
counted fractionally:

- per feature type and gene: `{mg,mt}.{ID,Name}.{CDS,rRNA,tRNA,ncRNA,tmRNA,CRISPR,…}_counts.tsv`;
- per functional term of each database: `{mg,mt}.mantis.{db}.counts.tsv` (e.g. `kegg_ko`, `kegg_module`, `cazy`).
  The terms are taken from the merged annotation, so they include all sources (Mantis consensus, Bakta), one line per
  gene and term (`scripts/split_annotation_terms.py`). Mantis reports Pfam domains by accession and by name
  (`PF00483` and `NTP_transferase` for the same domain); names are turned into accessions with the NAME / ACC pairs
  of `Pfam-A.hmm` (`annotation_terms: pfam_hmm`, default `<db>/mantis/pfam/Pfam-A.hmm`), earlier names of renamed
  families (e.g. `DUF541`, now `SIMPL`, PF04402) with Mantis' `pfam/metadata.tsv` and the DUF / UPF notes of
  `Pfam-A.hmm`, and each term is kept once per gene; names found in neither stay as they are; KEGG pathways given as `koNNNNN` and `mapNNNNN` are kept once as `koNNNNN`. Before 3.1.0 the counts came from
  Mantis' per-database GFFs, where a domain could appear under its accession and its name with different counts.

Coverage per gene and per contig (mosdepth, hts-nim-tools, BEDTools) goes to `Stats/{mg,mt}/`.

## Variants

`bcftools mpileup | bcftools call` on the metaG and metaT alignments: `Analysis/snps/{mg,mt}.variants.samtools.vcf.gz`.

## Protein database for metaproteomics

From the called proteins IMP3 builds a sample-specific search database for metaproteomics (length filters, cRAP
contaminants, optionally the proteomes of taxa detected by Kraken2): settings under `proteomics`.
