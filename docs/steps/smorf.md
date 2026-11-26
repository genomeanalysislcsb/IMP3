# smORF module

Small open reading frames (smORFs, ≤ 100 aa) are left out of the usual gene calling, annotation and metaproteomics
searches. The smORF module takes the small genes of the IMP3 annotation and predicts new ones in its intergenic
regions, removes false positives, counts the metagenomic and
metatranscriptomic reads on them, clusters them across samples and, optionally, builds per-sample metaproteomics
search databases with the IMP3 proteins, the smORFs, the human proteome and common contaminants.

Method and code are those of the smORF methods paper:
Davin, Ortís Sunyer, Delgado *et al.* High-resolution multi-omics enhances prediction and detection of
smORF-encoded proteins in the human gut microbiome. *Nat Commun* (2026),
[doi:10.1038/s41467-026-72762-5](https://doi.org/10.1038/s41467-026-72762-5); pipeline by Luis F. Delgado and
Júlia Ortís Sunyer
([smorf_methods_paper](https://gitlab.com/uniluxembourg/lcsb/systems-ecology/smorf_methods_paper)). Please cite
the paper when you use the module.

The module is an **independent analysis on top of finished IMP3 runs**: it reads the output folders of the samples
and writes to its own folder. It can run (again, e.g. with more samples or other settings) at any time after IMP3,
and it never changes the IMP3 results.

## Running it

In an [imp3 project](../portable.md), after `imp3 run`:

```bash
imp3 smorf -p <project> --dry-run      # what would run
imp3 smorf -p <project>                # local executor: on this machine; slurm: rules as Slurm jobs
imp3 smorf -p <project> --submit       # the Snakemake master as a Slurm job
imp3 smorf -p <project> --samples S01,S02 --outdir <project>/smorf_test
```

It uses every sample with `imp/<sample>/status/all.done` (or `--samples`) and writes to `<project>/smorf`. The first
run copies the settings to `<project>/smorf_config.yaml`; edit them there (for example `proteomics_db: build: true`
or `metadata:` for groups) and run again.

Without the launcher, for any folder of IMP3 output folders:

```bash
cp <IMP3>/config/config.smorf.yaml smorf.yaml      # set imp_root, outdir, db_path
snakemake -s <IMP3>/workflow/smorf/Snakefile --configfile smorf.yaml --use-conda --cores 24
```

**Software.** The rules use the conda environment `workflow/envs/IMP_smorf.yaml` (pyrodigal-gv 0.3.2, MMseqs2
15-6f452, Infernal 1.1.5, Aragorn 1.2.41, pyhmmer 0.10.14, subread, samtools; the versions of the paper). The
3.0 Apptainer image does not contain it yet: `imp3 smorf` creates it with conda in `conda_prefix` (or
`smorf_conda_prefix`, `<project>/conda_envs_smorf`), also in container projects. An image built with
`imp3 containerize` from an installation that has the environment contains it; set `smorf_container` in
`imp3.yaml` to use it.

**Databases** (in `<db_path>/smorf`, downloaded by the module when missing): Rfam covariance models (Rfam CURRENT,
~400 MB, `cmpress`ed), AntiFam 8.0, and for the proteomics databases the human RefSeq proteins (NCBI, deduplicated).
Other locations: `databases: rfam_cm / antifam_hmm / human_proteome`.

## What it does

Per sample, on the IMP3 assembly (`Assembly/<mgmt|mg|mt>.assembly.merged.fa(.gz)`) and its annotation
(`Analysis/annotation/annotation_CDS_RNA_hmms.gff(.gz)`). Where the smORF candidates come from
(`prediction: regions`):

- **`intergenic`** (default): the module does not call genes on the whole assembly again. The candidates are
  the CDS of the IMP3 annotation with ≤ 100 aa (`imp3_cds`; complete, not at a contig end) and new genes that
  pyrodigal-gv finds only in the intergenic regions of the IMP3 annotation: stretches of at least `min_region` nt
  without a feature of the types in `occupied` (CDS, tRNA, tmRNA, rRNA, ncRNA, CRISPR), extended by `flank` nt
  into the neighbouring features (default 0: new smORFs do not overlap IMP3 features). The genes found on these
  fragments are lifted back to contig coordinates. Since pyrodigal-gv sees only the fragments, it picks its
  metagenome model per fragment.
- **`full`**: the paper's pipeline, pyrodigal-gv on the whole assembly (e.g. to reproduce its numbers).

| step | tool | output (in the module's folder) |
|---|---|---|
| intergenic regions of the IMP3 annotation (`intergenic`) | | `Input/intergenic/<sample>.{fasta.gz,map.tsv}`, `Stats/<sample>.intergenic.tsv` |
| IMP3 CDS ≤ 100 aa (`intergenic`, `imp3_cds`) | | `pyrodigal_gv/imp3/<sample>.{smorf.faa,fna,gff}` |
| gene calling, run for genes ≥ 20 aa (`--min-gene 60 --max-overlap 60`) and for small genes (`--min-gene 15 --max-overlap 15`), metagenome mode, closed ends, no stop codon; on the intergenic regions or the whole assembly | pyrodigal-gv | `pyrodigal_gv/{default,small}/<sample>.{faa,fna,gff}` |
| smORF candidates of the gene calling: complete genes (`partial=00`) with an RBS motif, ≤ 100 aa (`prediction: max_aa`) | | `pyrodigal_gv/{default,small}/<sample>.smorf.faa` |
| identical candidates merged (100 % identity and coverage) | MMseqs2 | `potential_smORF/<sample>.smorf.{faa,fna,gff}` |
| ncRNA genes (Rfam, `--cut_tc`), tRNA / tmRNA genes (standard and bacterial code), spurious ORFs (AntiFam, gathering cut-offs), on all candidates | cmscan, Aragorn, pyhmmer | `potential_smORF/{RFAM,Aragorn_std,Aragorn_gcbact,antifam}/` |
| **final smORFs** (candidates without any of these hits) | | `smORF/<sample>.smorf.{faa,fna,gff}` |
| **IMP3 annotation with the smORFs** (`annotation_with_smorfs`): the IMP3 GFF, IMP3 CDS that are smORFs tagged `smORF=true`, new smORFs added as CDS (source `smORF_module`, ID `<sample>_smORF_<n>`, `product=small protein (predicted smORF)`) | | `smORF/<sample>.annotation_with_smORFs.gff` |
| numbers per step (intergenic bp, IMP3 vs new smORFs, filters, lengths) | | `Stats/<sample>.smorf_stats.tsv` |
| metaG / metaT reads per smORF on the IMP3 BAMs (`Assembly/<omic>.reads.sorted.bam`; unstranded, fractional multi-mapping / overlap as IMP3's gene counts); support class MG only / MT only / MG+MT / none | featureCounts | `read_support/<sample>.smorf_reads.tsv`, `read_support/<sample>.<omic>.counts.tsv` |

The IMP3 output itself is not changed: the annotation with the smORFs is a new file next to the smORFs, which
downstream tools (featureCounts, IMP explorer) can use in place of the IMP3 GFF.

smORF names are `<sample>::<set>_<contig>_<n>` (set `default`, `small` or `imp3`; `n` numbers the genes of a set
along the contig); `smORF/<sample>.smorf.gff` keeps the set prefix on the contig names, as in the paper's pipeline (`read_support/<sample>.smorf.saf` has the plain contig names).
`read_support: auto` counts reads for the samples whose BAMs are there; `true` requires them, `false` skips it.

All samples:

| step | output |
|---|---|
| clusters of all smORFs per level (`clustering: levels`, identity and coverage; default 100 % and 95 %) | `Cluster_predicted_smORF/Coverage_<c>_Ident_<i>_{rep_seq.fasta,cluster.tsv,rep_seq.fna}` |
| samples of each cluster, per group (`metadata`: TSV with sample and group) | `…_rep_samples.tsv`, `…_rep_presence_in_groups_report.tsv` |
| length distributions per group | `Size_distribution/smORF/*.png` |
| numbers: per sample (`smORF_samples.tsv`), per level: clusters, singletons, in one sample, in ≥ 2 samples, core (`smORF_clusters.tsv`), read support (`read_support_samples.tsv`) | |

## Metaproteomics search databases (`proteomics_db: build: true`)

As in the paper (MP_DB): per sample the Bakta proteins of the IMP3 run and the smORFs, deduplicated (100 % identity
and coverage), then together with the human RefSeq proteome and the CCP cRAP contaminants (14 June 2025, shipped in
`workflow/data/smorf/`), again deduplicated. Identical proteins in different samples get the same ID.

| output | content |
|---|---|
| `proteomics_db/final/<sample>.faa` | search database of the sample, headers `>source\|Cluster=<n>\|ID=<id>\|gene=…\|contig=…\|sample=…\|smORF=<True\|False>` |
| `proteomics_db/protein_clusters.tsv` | cluster, representative / member, protein ID, source (imp, smorf, human, crap), sample, contig, gene, ≤ 100 aa |

IDs are `<id_prefix>_<source>_<sample>_<n>` (`proteomics_db: id_prefix`, the paper used `expobiome`). Differences
to the paper's scripts: proteins are numbered per sample (the paper numbered all samples together), so samples
can be added later without renumbering; white space in headers is replaced by `_`.

## Licence

The scripts in `workflow/scripts/smorf/` come from the smORF methods paper and are under the GNU GPL v3
(`workflow/scripts/smorf/LICENSE.GPL-3.0`); the rest of IMP3 is under the MIT licence.
