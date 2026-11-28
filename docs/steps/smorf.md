# smORFs (after IMP3)

Small open reading frames (smORFs, ≤ 100 aa) are left out of the usual gene calling, annotation and metaproteomics
searches. `imp3 smorf` collects them for every finished sample of a project: the small genes of the IMP3 annotation
plus new small genes predicted in the parts of the assembly that carry no IMP3 feature. It runs any time after IMP3,
one job per sample (a single Python process, a few minutes per sample), and writes into the sample's IMP3 folder; the
IMP3 results themselves are not changed.

The gene calling follows the smORF methods paper: Davin, Ortís Sunyer, Delgado *et al.* High-resolution multi-omics
enhances prediction and detection of smORF-encoded proteins in the human gut microbiome. *Nat Commun* (2026),
[doi:10.1038/s41467-026-72762-5](https://doi.org/10.1038/s41467-026-72762-5); pipeline by Luis F. Delgado and
Júlia Ortís Sunyer ([smorf_methods_paper](https://gitlab.com/uniluxembourg/lcsb/systems-ecology/smorf_methods_paper)).
Please cite the paper when you use the module.

## Running it

In an [imp3 project](../portable.md), after `imp3 run`:

```bash
imp3 smorf -p <project> --dry-run           # the commands per sample
imp3 smorf -p <project>                     # every finished sample without status/smorf.done
imp3 smorf -p <project> --samples S01,S02   # only these
```

With `executor: slurm` every sample is one Slurm job (`smorf_threads`, default 16 threads; `smorf_mem_gb`,
`smorf_time`, `smorf_sbatch`), `parallel_samples` at a time (`-j`); with `executor: local` the samples run one after
the other. The first run copies the settings to `<project>/smorf_config.yaml` (edit them there) and downloads the
AntiFam HMMs to `<db_path>/smorf/AntiFam` (needs internet, e.g. the access node). Logs:
`imp_run_logs/smorf_<sample>_<job>.log`; `status/smorf.done` marks a finished sample (`--force` reruns it).

**Software.** `workflow/envs/IMP_smorf.yaml` (pyrodigal-gv 0.3.2, pyhmmer 0.10.14, subread, samtools; Infernal and
Aragorn for `paper: true`; the versions of the paper). `imp3 smorf` finds the environment in `conda_prefix` (or
`smorf_conda_prefix`, `<project>/conda_envs_smorf`) under the name Snakemake would give it, or uses `smorf_env` /
`smorf_container`; if it is missing, it prints the `micromamba create` command.

Without the launcher, for one IMP3 output folder:

```bash
python <IMP3>/workflow/scripts/smorf/imp3_smorf.py --sample-dir imp/S01 --config smorf_config.yaml \
    --db-path <db> --threads 16 [--outdir test_folder]
```

## What it does

| step | how |
|---|---|
| regions without a feature (`regions: intergenic`) | stretches of ≥ `min_region` nt without a feature of the types in `occupied` (CDS, tRNA, tmRNA, rRNA, ncRNA, CRISPR) of the IMP3 GFF, extended by `flank` nt; class `between_features`, `contig_end` (before the first / after the last feature) or `unannotated_contig` (contig without any feature) |
| IMP3 CDS ≤ `max_aa` aa (`imp3_cds`) | complete Bakta CDS not touching a contig end, translated with table 11 |
| new genes | pyrodigal-gv on the regions (`intergenic`) or the whole assembly (`full`), metagenome mode, closed ends, N runs masked; a run for genes ≥ 20 aa (`default`: min_gene 60, max_overlap 60) and one for small genes (`small`: 15, 15); genes ≤ `max_aa` aa, lifted to contig coordinates |
| exact duplicates | identical proteins merged (as MMseqs2 at 100 % identity and coverage in the paper); representative IMP3 > default > small, then by position; `copies` and `found_by` keep the rest |
| AntiFam | spurious ORFs (pyhmmer, gathering cut-offs) removed from the new candidates; Bakta already filters its own CDS |
| start status | `start_confirmed` / `start_open`, see below |
| RBS motif | reported, not filtered: `motif` / `none` / `not_assessable`; for IMP3 CDS from pyrodigal-gv on the gene plus 60 nt upstream (`unmatched` if it calls another start there) |
| overlaps | IMP3 features overlapping the smORF (other than the smORF itself), e.g. `regulatory_region:<riboswitch>` for leader peptides |
| metaproteomic status | `predicted` → `validated` with ≥ `min_unique_peptides` identified peptides that occur in the smORF and in no protein > 100 aa (IMP3 proteins, `metap_background`); I = L, modifications removed, any search engine |
| read support | featureCounts on `Assembly/<omic>.reads.sorted.bam` as IMP3's gene counts (`-O --fraction -M --fraction`, unstranded); class MG only / MT only / MG+MT / none |

`paper: true` with `regions: full` adds the paper's three further filters for a reproduction of its numbers: smORFs
without a Shine-Dalgarno motif, with an Rfam hit (`cmscan --cut_tc`, Rfam in `<db_path>/smorf/Rfam`) or an Aragorn
tRNA / tmRNA hit are removed. In intergenic mode Rfam and Aragorn would repeat what Bakta did on the whole assembly
(tRNAscan-SE, Aragorn, Infernal with Rfam); ncRNA genes and riboswitches are reported as overlaps instead, since
leader peptides (e.g. *trpL*) sit in riboswitch regions.

### Outputs

| file (in `imp/<sample>/`) | content |
|---|---|
| `Analysis/smORF/<sample>.smorf.faa`, `.fna`, `.gff` | the smORFs (IDs `<sample>_smORF_<n>`, contig coordinates) |
| `Analysis/smORF/<sample>.smorf_table.tsv` | one row per smORF: set, found_by, copies, position, length, region, start_status, rbs_status / motif / spacer, IMP3 locus tag, metaP status and peptides, overlaps, MG / MT reads, support |
| `Analysis/smORF/<sample>.smorf_stats.tsv` | numbers per step: regions per class, candidates per run, duplicates, filters, start and RBS status, new smORFs per region class, lengths, metaP, read support |
| `Analysis/smORF/<sample>.smorf.provenance.json` | settings, input files (size, time), database checksums (MD5), tool versions, run time |
| `Analysis/smORF/<sample>.smorf.<omic>.counts.tsv` | featureCounts tables |
| `Analysis/annotation/annotation_CDS_RNA_hmms_smORF.gff.gz` | the IMP3 GFF with the smORFs: IMP3 CDS that are smORFs get `smORF=true`, new smORFs are added as `CDS` (source `smORF_module`, `product=small protein (predicted smORF)`); attributes `smorf_set`, `smorf_region`, `smorf_start`, `smorf_rbs`, `smorf_status`, … (extra GFF3 attributes; other tools ignore them) |
| `status/smorf.done` | finished |

Across samples (clusters at 100 % and 95 % identity, presence per group, length distributions, read support and
metaP status per cluster) the smORFs are analysed by imp-explorer's cohort analysis.

## Complete genes, unconfirmed starts and RBS motifs

**All smORFs are complete genes in the formal sense.** pyrodigal-gv runs with closed ends (`-c`, as in the paper's
pipeline), and Bakta, the annotation of IMP3, also predicts genes with closed ends on linear contigs (every
metagenome contig; open ends only for circular replicons). A gene is therefore only called when its start codon and
its stop codon lie on the contig; every gene is `partial=00`, and the paper's "complete genes only" criterion
removes nothing.

**A start codon on the contig is not always the start of the gene.** ATG, GTG and TTG also occur inside genes.
When a contig begins in the middle of a gene, the real start lies outside the contig, and the gene caller reports
the gene from the next internal start codon instead:

```
real gene:      ATG.........................ATG..........TAA
contig:                   |-----------------ATG..........TAA------|
reported:                                   ATG..........TAA   (complete, partial=00)
```

The reported ORF is complete, but it is the 3' part of a longer gene; because it is short, it looks like a smORF.

**Start status.** The module looks upstream of the start codon in the same reading frame:

- `start_confirmed`: an in-frame stop codon lies upstream on the contig. A longer gene cannot run through it, so the
  gene begins at this start codon.
- `start_open`: no in-frame stop codon up to the contig end. The smORF may be a real small gene or the 3' part of a
  longer gene; this cannot be decided from the contig.

A stop codon occurs on average every ~21 codons in any reading frame of non-coding sequence or of another gene's
coding sequence, so genes with more than ~100 nt upstream on the contig are almost always `start_confirmed` — also
genes inside operons and the first gene of an operon (its upstream region, promoter and spacer, has stop codons in
all frames). `start_open` concerns genes near the 5' end of a contig, i.e. mostly short contigs. On BHC01 (median
contig 546 bp) 35 % of the paper pipeline's smORFs and 47 % of the IMP3 CDS ≤ 100 aa were `start_open`.

`start_open` smORFs are kept in the table and the GFF (`smorf_start=open`), but left out of the core set used for
statistics, clustering across samples (imp-explorer) and overlaps with sRNAs (`smorf_start=confirmed`). The same applies to all genes
of the IMP3 annotation near contig ends; the IMP3 gene counts do not use this status.

**RBS motif.** The paper keeps only smORFs for which Prodigal reports a Shine-Dalgarno motif (`rbs_motif` not
`None`). The module does not filter on it, but reports `rbs_motif`, `rbs_spacer` and whether at least `min_upstream`
(20) nt upstream were available (`motif`, `none`, `not_assessable`). The motif is missing for reasons other than a
false gene: genes near contig ends, and whole taxa that translate with weak or no Shine-Dalgarno motifs. All genes
of the pyrodigal-gv `default` run on the whole BHC01 assembly (434,536 genes; MAG phylum from GTDB-Tk of the IMP3
run):

| genes | without RBS motif |
|---|---|
| all | 34.8 % |
| > 100 aa / ≤ 100 aa | 31.7 % / 43.1 % |
| ≥ 20 nt / < 20 nt upstream on the contig | 29.6 % / 79.8 % |

| phylum (complete genes, ≥ 20 nt upstream) | genes | without motif | smORF-sized genes (≤ 100 aa) | without motif |
|---|---|---|---|---|
| Bacillota_A | 124,868 | 16.4 % | 15,308 | 20.0 % |
| Bacillota_I | 5,731 | 16.1 % | 938 | 20.9 % |
| Bacillota_C | 1,850 | 17.7 % | 226 | 13.7 % |
| Actinomycetota | 9,322 | 21.8 % | 1,077 | 27.5 % |
| Pseudomonadota | 5,945 | 26.1 % | 743 | 28.0 % |
| Bacteroidota | 32,373 | 36.6 % | 3,332 | 30.3 % |
| Desulfobacterota | 1,644 | 36.9 % | 240 | 32.5 % |
| Verrucomicrobiota | 4,133 | 81.9 % | 568 | 68.5 % |
| unbinned contigs | 203,457 | 36.4 % | 78,099 | 40.6 % |

An RBS filter would remove about 70 % of the Verrucomicrobiota (*Akkermansia*) smORFs, 30 % of the Bacteroidota
ones and 20 % of the Bacillota ones: a taxonomic bias, not a quality filter (large genes lack the motif as often in
these phyla). Apply it downstream (`rbs_status == motif`) if comparability with the paper is needed.

## Results on one Expobiome sample (BHC01)

Co-assembly of metaG + metaT, 463 Mb in 346,887 contigs (median 546 bp); IMP3 3.0.2, Bakta 1.9.x fork.

- **Port check.** The paper's original Snakemake pipeline and this module in `full` mode with the paper's filters
  gave identical results: 67,216 candidates and 66,757 final smORFs, same positions and protein sequences
  (comparison made with the former Snakemake version of the module, which shares the gene calling and filters).
- **Regions without a feature**: 600,580 regions, 162 Mb (35 % of the assembly): unannotated contigs 69.4 Mb
  (150,969, median 406 bp), contig ends 66.3 Mb (305,086, median 139 bp), between features 26.1 Mb (144,525,
  median 128 bp). Bakta calls genes with closed ends on linear contigs, so short contigs often carry none: CDS on
  29 % of the contigs < 500 bp, 64 % at 500–1,000 bp, 93 % at 1–2 kb and 100 % ≥ 2 kb.
- **Intergenic mode**: 119,457 smORFs = 112,403 IMP3 CDS ≤ 100 aa + 7,054 new (contig ends 4,755, between
  features 1,582, unannotated contigs 717); AntiFam removed 299; 8–100 aa, median 72 (former Snakemake version of the module, ~7 min on 24 cores).
- **Intergenic vs full (paper)**: 63,809 smORFs in both (96 % of the paper's set). Only in `full`: 2,948, of which
  2,257 lie inside or overlap an IMP3 CDS (excluded on purpose) and 690 in regions without a feature (missed by the
  fragment-wise gene calling). `flank: 60` recovers only 65 of them (flank 0: 19), but adds 5,647 new smORFs of which
  78 % overlap the end of an IMP3 CDS (and loses 1,367); `flank` therefore stays 0. Only in `intergenic`: 49,533 IMP3 CDS (no RBS filter on
  annotated genes) and 6,115 new genes.
- **Start status**: `start_open` 35 % of the paper's smORFs and 46 % of the intergenic ones (47 % of the IMP3 CDS
  ≤ 100 aa), mostly on short contigs.

