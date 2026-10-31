# Output files

Every sample has its own output folder (`outputdir`). The assembly prefix is `mgmt` for integrated runs (`mg` or
`mt` for single-omic runs). Sizes are those of one gut sample with ~35 M metagenome and ~47 M metatranscriptome
read pairs.

| folder | content | size |
|---|---|---|
| `Preprocessing/` | `{mg,mt}.{r1,r2,se}.preprocessed.fq.gz`; mt rRNA reads `mt.*.trimmed.rrna.fq.gz` | ~12 GB |
| `Assembly/` | `mgmt.assembly.merged.fa(.gz)` + BWA index, gene calls `mgmt.assembly.merged_1000.fa.frag.*`, `intermediary.tar.gz` | ~9 GB |
| `Analysis/annotation/` | `annotation_CDS_RNA_hmms.gff.gz`, Bakta (`bakta/`), Mantis, featureCounts tables | ~1.7 GB |
| `Analysis/taxonomy/` | `kraken/`, `bracken/`, `mOTUs/` | ~1.9 GB |
| `Analysis/snps/` | `{mg,mt}.variants.samtools.vcf.gz` | ~60 MB |
| `Binning/` | one folder per binner, `DASTool/`, `MetaWrap/` (final bins, `bins.contigs`), `per_bin_results.tar.gz` | ~0.8 GB |
| `Stats/` | read, mapping, coverage, bin statistics | ~1.3 GB |
| `Visualization/` | MetaWrap and Nonpareil plots | – |
| `status/` | `<step>.done` markers | – |

## Key tables

| file | content |
|---|---|
| `Stats/all_stats.tsv` | one table of the run: reads per preprocessing step, mapping, Kraken / mOTUs counts, annotation counts, assembly and bin numbers (`type`, `step`, `number`) |
| `Stats/all_bin_stats.tsv`, `Stats/mgmt/mgmt.bins.tsv` | final bins: CheckM2 completeness / contamination, N50, size, GC, GTDB-Tk classification and closest reference, GRiD |
| `Stats/{mg,mt}/{omic}.read_stats.tsv` | seqkit statistics of raw and trimmed reads |
| `Stats/{mg,mt}/{omic}.r{1,2}_nonpareil.diversity.txt` | Nonpareil sequence diversity |
| `Stats/{mg,mt}/mgmt.assembly.contig_read_count.txt` | reads per contig (contig, start, end, reads) |
| `Stats/{mg,mt}/mgmt.assembly.contig_depth.txt` | mean depth per contig |
| `Stats/{mg,mt}/mgmt.assembly.contig_flagstat.txt` | samtools flagstat of the reads on the assembly |
| `Stats/{mg,mt}/gene_read_count.txt`, `gene.mosdepth.*` | reads and coverage per gene |
| `Stats/mgmt/mgmt.assembly.length.txt`, `…gc_content.txt` | contig lengths and GC |
| `Analysis/annotation/{mg,mt}.mantis.{db}.counts.tsv` | featureCounts per functional term (`kegg_ko`, `kegg_module`, `kegg_pathway`, `cazy`, `enzyme_ec`, `pfam`, `go`, `cog`, `eggnog`, `tcdb`, `tigrfam`, …) |
| `Analysis/annotation/{mg,mt}.{ID,Name}.{CDS,rRNA,tRNA,ncRNA,CRISPR,…}_counts.tsv` | featureCounts per feature type, by ID or name |
| `Analysis/taxonomy/bracken/{omic}/{rank}.bracken` | Bracken estimates per rank (earlier runs: `kraken/{omic}.kraken_bracken_{rank}.report`) |
| `Analysis/taxonomy/mOTUs/{omic}.mOTU.{counts,rel_abun}.tsv` | mOTUs profiles |
| `Binning/MetaWrap/bins.contigs` | contig → final bin |

For a whole study, [IMP explorer](https://genomeanalysislcsb.github.io/imp-explorer/) collects these tables from all
samples into one release and an interactive app.
