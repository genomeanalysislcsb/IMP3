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

Kraken2 settings: `kraken.krakendb`, `kraken.confidence` (default 0.501), `kraken.min_hit_groups` (3). The database
decides the names (NCBI or GTDB): see [Databases](../databases.md#choosing-the-kraken2-database). An experimental
dual-database mode classifies against a GTDB and an NCBI nt database in one run.

By default Kraken2 loads its database into memory: the threads are sized from the database size and the
big-memory specs (`cluster_specs`), and on ULHPC the rules go to iris bigmem (`profiles/slurm/cluster.ulhpc.yaml`).
With `kraken.memory_mapping: true` it runs with `--memory-mapping` (the database is read through the page cache)
on a whole regular node instead. That only pays off when the database fits into a regular node's memory, where it
stays cached for the following Kraken2 jobs; a larger one (core_nt, ~280 GB, vs ~210 GB per aion node) is re-read
from disk again and again. On ULHPC, also remove the iris bigmem entries of the Kraken2 rules from the project's
`profile/cluster.yaml`.
