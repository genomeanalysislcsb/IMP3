# Many samples, any machine: `imp3`

`imp3` (in the IMP3 folder) runs IMP3 for every sample of a project. The same command works in four setups:

| | conda environments on the host | IMP3 container (Apptainer / Singularity) |
|---|---|---|
| **Slurm cluster** (every rule a job, bigmem and GPU rules on their partitions) | `executor: slurm` | `executor: slurm` + `container:` |
| **One machine** (server or workstation, no scheduler) | `executor: local` | `executor: local` + `container:` |

It replaces the per-site submit scripts (`IMP_parallel_submit.sh` and its templates): it finds the reads of every
sample, writes the sample's IMP3 config from a template, skips finished samples (`status/all.done`), starts one
Snakemake master per sample (several at a time) and keeps a numbered log per run.

Requirements on the host: the Snakemake 7.32 environment of IMP3 (`requirements.yaml`, with the
[shadow patch](installation.md#snakemake-7-shadow-patch)); for the container also `apptainer` or `singularity`.
The databases stay outside the container in every setup.

## A project

```bash
/path/to/IMP3/imp3 init my_project --executor slurm --site ulhpc --db /path/to/imp3_db \
    --conda-prefix /path/to/imp3_envs --config-template /path/to/IMP_submit_template_config.yaml
```

```text
my_project/
├── imp3.yaml               settings (below)
├── sample_config.yaml      IMP3 config template (__SAMPLE__, __IMP_OUT__, __IMP_DB__, __MG_READS__, __MT_READS__)
├── profile/                Snakemake profile: config.yaml (+ cluster.yaml for Slurm)
├── samples/<sample>/       read files of each sample (or symlinks); names with mg / mt, R1 before R2
├── config/                 generated config per sample (imp_config_<sample>.yaml)
├── imp/<sample>/           IMP3 output
├── imp_run_logs/           imp_<sample>_001.log, _002.log, ... (one per run)
└── slurm_logs/             logs of masters submitted with --submit
```

Without `--config-template`, the project gets `profiles/sample_config.template.yaml` (all settings of the current
IMP3 version; Mantis and proteomics databases below `__IMP_DB__`). Instead of `samples/`, a sample sheet
(`sample_sheet:` TSV with `sample`, `mg_r1`, `mg_r2`, `mt_r1`, `mt_r2`) can list the reads.

```bash
cd my_project
/path/to/IMP3/imp3 run --dry-run            # configs + Snakemake dry run per sample
/path/to/IMP3/imp3 run                      # all unfinished samples, parallel_samples at a time
/path/to/IMP3/imp3 run --samples GF06,GF12  # only these
/path/to/IMP3/imp3 run -- --until assembly  # extra Snakemake options after --
/path/to/IMP3/imp3 status                   # done / running / failed / stopped per sample, steps done, last log
```

On a cluster, `imp3 run` runs the Snakemake masters on the machine where it is started (e.g. in tmux on a login
node, as before). `imp3 run --submit` instead submits every master as a long Slurm job (`submit_master` options,
1 core), so nothing has to keep running on the login node; each master then submits the rules as usual.

### Settings (`imp3.yaml`)

| key | default | meaning |
|---|---|---|
| `imp3_dir` | the IMP3 folder of `imp3` | IMP3 checkout (Snakefile) |
| `db_path` | – | IMP3 databases (`__IMP_DB__`) |
| `reads` / `sample_sheet` | `samples` / – | where the reads are |
| `config_template` | `sample_config.yaml` | IMP3 config template |
| `executor` | `slurm` | `slurm` or `local` |
| `profile` | `profile` | Snakemake profile |
| `conda_prefix` | – | conda environments of IMP3 (native runs) |
| `container` | – | IMP3 `.sif`: every rule runs in the container |
| `bind` | – | extra folders the container must see (reads, databases, IMP3 and the project are bound automatically) |
| `container_mounts` | – | mount points built into the image (see *Mount points*); folders below them are bound as the mount point |
| `parallel_samples` | 10 (local: 1) | Snakemake masters at the same time |
| `submit_master` | `--time 5-00:00:00 --qos long --cpus-per-task 1 --mem 8G` | `sbatch` options for `imp3 run --submit` |
| `pre_commands` | – | shell lines before Snakemake, e.g. `module load tools/Singularity` |
| `tmp_root` | `/dev/shm` (local: `/tmp`) | root of IMP3's temporary folder |
| `cores`, `mem_gb`, `gpu` | all, 90 %, auto | local: what IMP3 may use on the machine |
| `snakemake_args` | – | extra Snakemake options for every run |

## On a Slurm cluster

`imp3 init --executor slurm --site <site>` writes `profile/config.yaml` (cores and jobs per sample, `status.py` for
the job state, `use-conda`) and `profile/cluster.yaml`, the `sbatch` options per rule:

- `--site ulhpc`: aion for most rules, iris `bigmem` for Kraken2, an iris volta GPU for COMEBin training
  (Slurm multi-cluster);
- `--site generic`: one cluster with the partitions `batch`, `bigmem` and `gpu` (`--gres gpu:1`) and empty `qos` /
  `account`: adapt them to your cluster.

```yaml
__default__:
  call: "sbatch"
  partition: "--partition batch"
  qos: ""
  account: "--account myproject"
  mem_per_cpu: "--mem-per-cpu "     # IMP3 asks memory per core (cluster_specs of the sample config)
  ...
kraken_reads:
  partition: "--partition bigmem"
comebin_train:
  partition: "--partition gpu"
  threads: "--gres gpu:1 --cpus-per-task "
```

Memory: every rule asks `mem_mb` per core from `cluster_specs` of the sample config (`main_cluster` for normal
nodes, `big_mem_*` for large-memory rules, `high_mem_cluster` for a separate large-memory cluster). Set them to the
nodes of your cluster.

## On one machine

```bash
/path/to/IMP3/imp3 init my_project --executor local --db /path/to/imp3_db --conda-prefix /path/to/imp3_envs
```

`imp3 run` then starts Snakemake with `--cores` (all, or `cores:`) and replaces `cluster_specs` by the machine: every
rule gets the same memory per core, so cores × memory per core never exceeds the machine (`mem_gb`, default 90 % of
the RAM). With a GPU (`nvidia-smi`; `gpu: 0` to switch off) COMEBin trains on it (`--resources gpu=N`, one COMEBin
at a time per GPU); without one COMEBin runs on the CPU (slow), or drop it from `binning: binners`.

Samples run one after the other, each with all cores. The largest Kraken2 databases need ~450 GB of memory: on a
smaller machine use a smaller Kraken2 database.

## The IMP3 container

One Apptainer/Singularity image holds every conda environment of IMP3, as installed and tested on the cluster
(including the patched COMEBin, Bakta and quickmerge and the mOTUs database), at `/conda-envs/<hash>`, where
Snakemake looks for them in a containerized image. Snakemake itself runs on the host; every rule runs in the
container (`--use-conda --use-singularity`). Nothing is installed at run time.

```yaml
# imp3.yaml
container: /path/to/imp3.sif
pre_commands: ["module load tools/Singularity"]     # if apptainer / singularity is a module
```

`imp3 run` then

- sets `container:` in the sample config (the Snakefile declares it as `containerized:`);
- binds the IMP3 folder, the project, the databases, the read folders and every existing path of the config
  (Mantis folders, proteomics files), plus `bind:`;
- keeps these paths as written (only symlinked files are followed, through every link of a chain, each folder bound): on a site with several clusters, use paths that
  exist on all of them, e.g. `/work/projects/...` on ULHPC rather than `/mnt/aiongpfs/...` (aion only);
- binds each as `path:path`: without a destination Apptainer resolves symlinks in the path and mounts the folder
  at its real path, so e.g. `/work/projects/...` would not exist in the container;
- checks once, with the image, that every bound folder is visible in the container, and stops a sample that would
  not find its files (see *Mount points* below); without apptainer on the node it only says so;
- checks that the temporary folder (`tmp_root`, default `/dev/shm`: IMP3's temporary files and Snakemake's shadow
  folders) is the same inside the container. Some sites give the container its own empty `/dev/shm` (ULHPC), and
  shadow rules then fail with empty logs; set `tmp_root` to a shared file system there, e.g.
  `/scratch/users/$USER/imp3_tmp` (`imp3 init --site ulhpc --container ...` does);
- adds `--nv`, so COMEBin sees the GPU of a GPU node (a warning on nodes without GPU, harmless);
- creates the marker files of IMP3's setup rules (`submodules/comebin.ready`, ...) if missing: the setup was done when
  the image was built.

On Slurm the jobs are submitted as before and wrap their command in `apptainer exec`; the compute nodes need
apptainer / singularity.

### Building the image

The image is built from the conda environments of an existing installation, e.g. the cluster where IMP3 runs, so it
contains exactly the tested environments:

```bash
/path/to/IMP3/imp3 containerize --envs /path/to/imp3_envs -o imp3_container
cd imp3_container && ./build.sh
```

`imp3 containerize` finds the environment of every `workflow/envs/*.yaml` (and of the binny module in
`submodules/binny`) under the prefix (Snakemake's hash of the prefix and the environment file), lists them in
`envs.tsv` and writes `pack_envs.sh`, `imp3.def` (miniforge base, the environments, tools for the rules without an
environment) and `build.sh`. The environments are copied unchanged (tar) to their original path in the image, e.g.
`/mnt/.../envs/<hash>`, and linked from `/conda-envs/<hash>`: nothing inside them has to be rewritten, so they are
exactly the tested ones (relocating them with conda-pack needs the original package caches, which older
installations often no longer have). At run time, do not bind a host folder over that original path at another
site; on the cluster itself the same environments are there, so it does no harm. Environments that are not
installed (e.g. binners never used) are reported and left out; a rule that needs one fails in the container.

`build.sh` needs apptainer with `--fakeroot` (or root; `BUILD_ARGS=` to change), and free space of about twice the
environments in the build folder and `APPTAINER_TMPDIR`. The image is large (tens of GB); build it once and copy
the `.sif` to other sites. It must belong to the same IMP3 version: the environment hashes come from the files in
`workflow/envs`.

`build.sh` also writes `imp3_submodules.tar.gz`: the third-party code in `submodules/` that Snakemake reads on the
host (the binny module and its marker database, the Bakta and COMEBin forks, ESOM). It is not part of the
repository; at another site unpack it into the IMP3 checkout of the same version
(`tar -xzf imp3_submodules.tar.gz -C /path/to/IMP3`), next to the `.sif`.

### Mount points

Many Apptainer installations can create a bind destination that does not exist in the image (overlay); others,
e.g. on ULHPC, cannot and skip such a bind without an error, so the rules do not find their files. There, the
image needs the top folders of the site's file systems, and IMP3 binds those instead of the single folders:

```bash
imp3 containerize --envs ... --mounts /work/projects,/scratch/users     # created in the image (label imp3.mounts)
```

```yaml
# imp3.yaml
container_mounts: [/work/projects, /scratch/users]     # imp3 init --site ulhpc --container ... sets these
```

A folder below a mount point is then bound as the mount point (e.g. `/work/projects`); other folders are bound as
they are. The mount points must not contain the environments' original path (`containerize` refuses that). An
image built at one site carries that site's mount points; at another site without overlay, add yours in a thin
layer on top of it (a few minutes, the environments are not repacked):

```bash
cat > mounts.def <<'EOF'
Bootstrap: localimage
From: imp3.sif

%post
    mkdir -p /data/projects /scratch
EOF
apptainer build --fakeroot imp3_site.sif mounts.def
```

and set `container` and `container_mounts` in `imp3.yaml` accordingly. `apptainer inspect --labels imp3.sif` shows
the mount points of an image. Sites with overlay need none of this.

## Notes

- IMP3's temporary folder is `/dev/shm/<run>` by default (RAM). On machines with little shared memory set
  `tmp_root` (`imp3.yaml` or the sample config), e.g. `/tmp` or a local scratch.
- `imp3 run` always unlocks a sample before starting it (as the old submit scripts did): do not start the same
  sample twice.
- The old scripts in `runscripts/` still work; `imp3` does not change how IMP3 itself runs.
