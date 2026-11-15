# FAQ and troubleshooting

**A conda environment does not resolve.**
Set `channel_priority: flexible` and use only conda-forge and bioconda in `~/.condarc`
([Installation](installation.md#2-snakemake-environment)).

**Outputs of a shadow rule are missing after the job finished ("Unshadowing failed").**
Snakemake 7 needs the [shadow patch](installation.md#snakemake-7-shadow-patch).

**Kraken2, GTDB-Tk or the assembly run out of memory.**
Large Kraken2 databases need up to ~450 GB; GTDB-Tk ~100 GB. Set the big-memory nodes in `cluster_specs`, choose a
smaller Kraken2 database, or run these jobs on a large-memory partition (cluster config).

**Mantis fails with "oom_kill" although its outputs were written.**
Slurm reports the job as `OUT_OF_MEMORY` (`Detected … oom_kill events` at the end of
`logs_cluster_jobs/IMP3.mantis_annot.<jobid>.stdout`), Snakemake counts it as failed, removes the outputs and retries.
Up to IMP3 v3.0.1 Mantis worked in `/dev/shm`, which is memory (tmpfs) and counts against the job's limit together
with Mantis itself; the intermediate and compressed files pushed the job over it, on large-memory nodes too. Since
v3.0.2 Mantis works in a job-private folder in the shadow directory (below `tmp_root`, on disk): update the IMP3
checkout and restart the sample (`imp3 run --samples <sample> --submit`); finished steps are kept.
If Mantis still runs out of memory on very large protein sets, run `mantis_annot` on large-memory nodes in
`profile/cluster.yaml` (ULHPC: uncomment the block below "Only if aion oom"; other clusters: your big-memory partition)

```yaml
mantis_annot:
  clusters: "--clusters 'iris'"        # ULHPC only
  partition: "--partition bigmem"
```

The annotation does not depend on the node.

**The run stopped; how do I resume?**
Run the same command again with `--rerun-incomplete`. Finished steps (`status/*.done`) are not repeated.

**A Snakemake lock blocks the run** ("Directory cannot be locked").
Make sure no other Snakemake runs on this output folder, then `snakemake … --unlock` once.

**Which Kraken2 database names do I get?**
The names of the database you installed: NCBI for standard / PlusPFP / nt, GTDB for GTDB-based databases
([Databases](databases.md#choosing-the-kraken2-database)).

**Can I run only part of the workflow?**
Yes: set `steps` and give the matching [inputs](inputs.md), e.g. only `analysis binning taxonomy summary` on an
existing assembly.

**Where do I ask questions or report bugs?**
Open an issue on [GitHub](https://github.com/genomeanalysislcsb/IMP3/issues).
