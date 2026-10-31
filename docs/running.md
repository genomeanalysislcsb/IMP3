# Running IMP3

## Locally or in one job

```bash
conda activate snakemake_env
snakemake --snakefile Snakefile --configfile config/my_sample.yaml \
          --use-conda --conda-prefix /path/to/conda_envs --cores 28 --reason --printshellcmds
```

Always start with a dry run (`--dry-run`).

## On a Slurm cluster (ULHPC aion / iris)

The profile in `cluster_profile/config.yaml` submits every rule as its own job:

```bash
snakemake --profile cluster_profile/ --configfile config/my_sample.yaml
```

Set in `cluster_profile/config.yaml` (lines marked `### USER_INPUT ####`): `conda-prefix`, `cluster-config`
(absolute path to `config/slurm_simple_ul_hybrid_cluster.config.yaml`), `cluster-status`
(`runscripts/status.py`), cores and jobs.

## Many samples

`runscripts/` holds helpers used on the ULHPC:

| script | purpose |
|---|---|
| `runIMP3_submit.sh` | submit one IMP3 run per sample of a sample list |
| `runIMP3_tmux.sh` | run samples from a tmux session on a login node |
| `launch.imp.sh`, `IMP3_withReport.sh`, `IMP3_withoutReport.sh` | single-sample launchers |
| `run_bracken.sh` | rerun Bracken on existing Kraken2 reports |
| `run_manual_imp3_mapping.sh`, `run_manual_imp3_featurecounts.sh` | redo the mapping / featureCounts of a sample |
| `status.py` | Slurm job status for Snakemake |

A sample list is a TSV with `Sample`, `MG_R1`, `MG_R2`, `MT_R1`, `MT_R2`.

## Re-running and resuming

Snakemake reruns only what is missing or outdated; `status/*.done` marks finished steps. After an interrupted run,
add `--rerun-incomplete`. To redo a step, delete its `status/<step>.done` and the outputs, or use
`--forcerun <rule>`.
