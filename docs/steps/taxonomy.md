# Taxonomy

The community composition is profiled from the reads (Kraken2 + Bracken, mOTUs), the contigs are classified, and the
MAGs get a GTDB taxonomy (see [Binning](binning.md)). Output: `Analysis/taxonomy/`.

| method | input | output |
|---|---|---|
| [Kraken2](https://github.com/DerrickWood/kraken2) | preprocessed reads of each omic (pairs and singletons; reports combined) | `kraken/{omic}.kraken.report`, `…unfiltered.report` |
| [Bracken](https://github.com/jenniferlu717/Bracken) | Kraken2 reports, per rank (`bracken.levels`, default domain … species) | `bracken/{omic}/{rank}.bracken` |
| Kraken2 on contigs | the assembly; per-contig taxonomic entropy | `kraken/mgmt.kraken.contig.*` |
| [mOTUs](https://github.com/motu-tool/mOTUs) | preprocessed reads | `mOTUs/{omic}.mOTU.{counts,rel_abun}.tsv` |
| mOTU links | marker genes on the contigs (fetchMG) linked to mOTUs | `mOTU_links/` |
| [EukDetect](https://github.com/allind/EukDetect) (optional) | reads against BUSCO genes of eukaryotes | `eukdetect/` |

Kraken2 settings: `kraken.krakendb`, `kraken.confidence` (default 0.501), `kraken.min_hit_groups` (3). The database
decides the names (NCBI or GTDB): see [Databases](../databases.md#choosing-the-kraken2-database). An experimental
dual-database mode classifies against a GTDB and an NCBI nt database in one run.
