# Configuration

A run is defined by one YAML file; `config/config.imp.yaml` is the template (`config/config_test.yaml` for the
test sample). Lines marked `### USER_INPUT ####` must be set for every installation or project. The schema is in
`schemas/schema.imp.yaml`.

- Only `raws` and `outputdir` are required; every other field has a default.
- The order of the fields does not matter.
- IMP3 writes the complete configuration it used, defaults included, to `<outputdir>/sample.config.yaml`.

## Run

| key | example | meaning |
|---|---|---|
| `steps` | `"preprocessing assembly analysis binning taxonomy summary"` | steps to run |
| `sample` | `sample01` | sample name, prefixed to the contig names; punctuation becomes `_`; empty = the last two parts of `outputdir` |
| `outputdir` | `/path/to/IMP3_outputs/sample01` | output folder of the sample |
| `tmp_dir` | `/path/to/tmp` | temporary files |
| `db_path` | `/path/to/IMP3_databases` | databases from the install workflow |
| `raws.*` | | input files, see [Inputs](inputs.md) |
| `summary_steps` | `"stats"` | `stats`, `vis` |
| `compress_level` | `mid` | how much of the intermediate output is compressed at the end |

## Resources

`cluster_specs` describes the machines: memory per core and in total, cores, for normal and big-memory nodes
(`main_cluster`, `high_mem_cluster`); IMP3 derives threads and memory per rule from it. `mem_unit` is `K`, `M`, `G`
or `T`.

## Preprocessing

| key | default | meaning |
|---|---|---|
| `fastp.adapter_trimming.{mg,mt}` | auto-detection | adapter sequences or FASTA |
| `fastp.quality_filtering` | Q15, ≤ 40 % low-quality bases, ≤ 5 N | read filtering |
| `fastp.global_trimming` | front 10, tail 1 | fixed trimming of R1 (R2 follows R1 unless set) |
| `fastp.quality_cutting` | tail and right, window 1, Q20 | sliding-window trimming |
| `fastp.length_filtering.length_required` | 40 | minimum read length |
| `fastp.merging.merge` | true | merge overlapping pairs |
| `preprocessing_filtering`, `filtering.filter` | `PhiX T2TCHM13v2genomic T2TCHM13v2transcripts` | references whose reads are removed |
| `keep_rrnas` | false | keep rRNA reads in the preprocessed mt reads |
| `sortmerna.files` | SILVA / Rfam rRNA sets | rRNA databases |
| `read_coverage_normalization` | off | normalise mt read coverage |

## Assembly

| key | default | meaning |
|---|---|---|
| `assembly.hybrid` | true | assemble mg and mt reads together |
| `assembly.assembler` | `megahit` | `megahit` or `metaspades` |
| `assembly.merge` | `quickmerge` | merging of the assembly rounds |
| `assembly.kmer_list`, `mink`, `maxk`, `step` | 21–141 | k-mers |

## Annotation and functions

| key | default | meaning |
|---|---|---|
| `annotation` | `mantis` | functional annotation (gene calls by Bakta) |
| `mantis.default.*_ref_folder` | KOfam | reference folders (KOfam, NCBI, eggNOG, Pfam, TCDB) |
| `hmm_DBs`, `hmm_settings` | | extra HMM databases and their cut-offs |
| `featureCountsStranding` | `mt: 2`, `mg: 0` | strandedness for featureCounts |

## Taxonomy

| key | default | meaning |
|---|---|---|
| `kraken.krakendb` | e.g. `k2_nt_20231129` | Kraken2 database (NCBI nt, standard, or GTDB-based) |
| `kraken.confidence` | 0.501 | Kraken2 confidence |
| `kraken.min_hit_groups` | 3 | minimum hit groups |
| `bracken.levels` | `D P C O F G S` | Bracken ranks |
| `motus`, `eukdetect` | | mOTUs and EukDetect settings |

## Binning

| key | default | meaning |
|---|---|---|
| `binning.binners` | e.g. `MetaDecoder` | binners to run (CONCOCT, MetaBAT, MetaDecoder, SemiBin, COMEBin, binny, MaxBin) |
| `binning.<binner>.cutoff` | 1000–2000 | minimum contig length per binner |
| `binning.DAS_Tool.score_threshhold` | 0.6 | DAS Tool score |
| `binning.MetaWrap.min_completeness`, `max_contamination` | 60, 10 | MetaWrap refinement |
| `binning.binny.*` | | binny settings |

## Proteomics (optional)

`proteomics` builds a sample-specific protein database for metaproteomics from the predicted genes (lengths,
cRAP contaminants, proteomes of Kraken2-detected taxa).
