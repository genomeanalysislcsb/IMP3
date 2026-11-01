# Running IMP3

!!! tip "Many samples"
    For many samples, on a Slurm cluster or one machine, with or without the IMP3 container, use the
    [`imp3` command](portable.md). This page shows the underlying Snakemake calls.

## Locally or in one job

```bash
snakemake -s /path/to/IMP3/Snakefile --configfile sample.yaml \
          --use-conda --conda-prefix /path/to/IMP3/conda --cores 28
```

| option | meaning |
|---|---|
| `--configfile` | the [config](configuration.md) of the sample |
| `--use-conda --conda-prefix` | conda environments of IMP3; use the same prefix as at installation |
| `--cores` | maximum cores used at the same time (memory per rule comes from `cluster_specs`) |
| `--dry-run` | show what would run |
| `--report report.html` | after a run: report with environments, rules and run times |

Always start with a dry run.

## On a Slurm cluster

The profile in `cluster_profile/config.yaml` submits every rule as its own job:

```bash
snakemake --profile cluster_profile/ --configfile config/my_sample.yaml
```

Set in `cluster_profile/config.yaml` (lines marked `### USER_INPUT ####`): `conda-prefix`, `cluster-config`
(absolute path to `config/slurm_simple_ul_hybrid_cluster.config.yaml`), `cluster-status`
(`runscripts/status.py`), cores and jobs.

## Many samples

`runscripts/` holds helper scripts for Slurm clusters (adapt the paths inside them):

| script | purpose |
|---|---|
| `runIMP3_submit.sh` | submit one IMP3 run per sample of a sample list |
| `runIMP3_tmux.sh` | run samples from a tmux session on a login node |
| `launch.imp.sh`, `IMP3_withReport.sh`, `IMP3_withoutReport.sh` | single-sample launchers |
| `run_bracken.sh` | rerun Bracken on existing Kraken2 reports |
| `run_manual_imp3_mapping.sh`, `run_manual_imp3_featurecounts.sh` | redo the mapping / featureCounts of a sample |
| `status.py` | Slurm job status for Snakemake |

A sample list is a TSV with `Sample`, `MG_R1`, `MG_R2`, `MT_R1`, `MT_R2`.

## Other batch systems

Without the profile, give the submission command and a cluster config file directly:

```bash
snakemake -s Snakefile --configfile sample.yaml --use-conda --conda-prefix conda --cores 28 \
  --cluster-config config/slurm.config.yaml \
  --cluster "{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem_mb} {cluster.threads}{threads} {cluster.partition} {cluster.nodes}"
```

The cluster config (YAML) holds the submission command and the argument names, with a default section and
sections for rules that need large-memory nodes. Templates: `config/slurm.config.yaml`,
`config/slurm_simple.config.yaml`, `config/slurm_simple_ul_hybrid_cluster.config.yaml` (example for a cluster with normal and large-memory nodes),
`config/pbs.config.yaml` and `config/sge.config.yaml`.

```yaml
__default__:
  call: "sbatch"
  nodes: ""
  mem_per_cpu: "--mem-per-cpu "
  partition: "-p batch"
  runtime: "-t"
  threads: "-n"
```

## Re-running and resuming

Snakemake reruns only what is missing or outdated; `status/*.done` marks finished steps. After an interrupted run,
add `--rerun-incomplete`. To redo a step, delete its `status/<step>.done` and the outputs, or use
`--forcerun <rule>`.
