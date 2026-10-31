# Summary

The summary step gathers the statistics of all steps into two tables and archives intermediate files.

| file | content |
|---|---|
| `Stats/all_stats.tsv` | reads per preprocessing step, mapping, Kraken2 / mOTUs counts, annotation counts, assembly statistics, number and quality of bins (`type`, `step`, `number`) |
| `Stats/all_stats.Rdata` | the same as an R workspace |
| `Stats/all_bin_stats.tsv` | one row per final bin |

`summary_steps: "stats"` is the default; the plotting step (`vis`) is currently disabled. `compress_level`
decides how much intermediate output is compressed or removed at the end.

For a whole study, the [IMP explorer](../downstream.md) collects these tables (and the taxonomy, bin abundance and
function counts) from all samples into one interactive app.
