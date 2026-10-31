# Databases

All databases live in one folder, `db_path`, set in both `config/config.install.yaml` and the run config. The install
workflow (`install/Snakefile`) downloads and prepares them; some are downloaded on first use by the main workflow.

| used by | database | installed by | config key |
|---|---|---|---|
| preprocessing | PhiX genome (and host genomes, e.g. human T2T-CHM13 genome + transcripts) | `genomes` in `config.install.yaml` | `filtering.filter` |
| preprocessing (metaT) | SortMeRNA rRNA sets (SILVA 16S/18S/23S/28S, Rfam 5S/5.8S) | install workflow | `sortmerna.files` |
| analysis | Bakta database | main workflow on first use (`get_bakta_db`) | – |
| analysis | Mantis references (KOfam, Pfam, eggNOG/NOG, NCBI, TCDB) | install workflow (`mantis_setup`) | `mantis.default.*_ref_folder` |
| analysis | extra HMM collections | by hand | `hmm_DBs`, `hmm_settings` |
| taxonomy | Kraken2 database (e.g. MiniKraken, PlusPFP, NCBI nt, GTDB) + Bracken k-mer files | install workflow (`kraken2`) | `kraken.krakendb` |
| taxonomy | mOTUs database | install workflow / first use | `motus` |
| taxonomy | EukDetect BUSCO database (optional) | by hand | `eukdetect.database_dir` |
| binning | CheckM2 database | main workflow on first use (`get_checkm2_db`) | – |
| binning / taxonomy | GTDB-Tk reference data | install workflow (`gtdbtk`) | – |

## Space and memory (approximate)

| database | disk | memory to use it |
|---|---|---|
| MiniKraken2 8 GB / PlusPFP-8 | 8 GB | ~10 GB |
| Kraken2 standard / NCBI nt | 50 GB – 500 GB | up to ~450 GB (or use memory mapping) |
| GTDB-Tk | ~100 GB | ~100 GB (pplacer) |
| Mantis (all references) | ~100–200 GB | moderate |

Choose small databases for a first test; the test run works with MiniKraken2.

## Choosing the Kraken2 database

The Kraken2 database decides the names in the read-based taxonomy: NCBI-based databases (standard, PlusPFP, nt)
give NCBI names, GTDB-based databases give GTDB names, which then match the GTDB-Tk taxonomy of the MAGs. IMP3 also
has an experimental dual-database mode that runs reads and contigs against a GTDB and an NCBI nt database
(`run_dual_db`).
