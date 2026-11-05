# Installation and test run

IMP3 needs Linux, conda (or mamba / micromamba), Snakemake 7 and, depending on the databases, up to ~1 TB of disk
space and ~450 GB of memory for the largest Kraken2 databases. All tools run in conda environments that Snakemake
creates; third-party code that is not on conda goes to `submodules/` (not part of the repository): see
[Third-party code](#third-party-code-submodules).

Run all commands from the IMP3 folder.

## 1. Clone the repository

```bash
git clone https://github.com/genomeanalysislcsb/IMP3.git
cd IMP3
```

## 2. Snakemake environment

```bash
conda env create -n snakemake_env -f requirements.yaml     # Snakemake 7.32, pigz, filetype
conda activate snakemake_env
```

Or add the two extra packages to an existing Snakemake 7 environment:
`conda install pigz=2.6 filetype=1.2.0 -c conda-forge`.

Recommended `~/.condarc`, so that all environments resolve:

```yaml
channel_priority: flexible
channels:
  - conda-forge
  - bioconda
```

!!! warning "Snakemake 7 and the shadow directive"
    IMP3 uses Snakemake's `shadow` directive. With Snakemake 7 the function `unshadow_output` in Snakemake's
    `dag.py` must be patched so that it copies nested outputs back correctly; see
    [the patch](#snakemake-7-shadow-patch) at the end of this page.

## 3. Databases and tools

Edit `config/config.install.yaml`:

- `db_path`: where the databases go (up to ~1 TB, depending on the Kraken2 and Mantis databases);
- `kraken2`: the Kraken2 database(s) to download (comment out the others);
- `mantis`: the Mantis reference databases;
- `binny` › `snakemake_env_path`: path of your Snakemake environment (`echo $CONDA_PREFIX` while it is active), or
  `in_path` if Snakemake is always on the `PATH`.

Then run the install workflow (dry run first). The conda environments go to `--conda-prefix`; use the same prefix
for all later runs.

```bash
conda activate snakemake_env
snakemake -s install/Snakefile --configfile config/config.install.yaml \
          --use-conda --conda-prefix ./conda --cores 1 --dry-run
snakemake -s install/Snakefile --configfile config/config.install.yaml \
          --use-conda --conda-prefix ./conda --cores 14 --printshellcmds
```

This may take hours for large databases. What is installed is listed under [Databases](databases.md).

## 4. Test run

The test reads are in `test/Reads/` (metaG and metaT, paired, gzipped). In `config/config_test.yaml` set:

- `cluster_specs`: memory and cores of your machine or cluster nodes (set the `big_mem` values equal to the normal
  ones if you have no large-memory nodes);
- `tmp_dir`, `outputdir`, `db_path`, `kraken.krakendb` and the Mantis reference folders, as installed in step 3;
- `raws`: the test reads, `test/Reads/MG.R{1,2}.fq.gz` and `test/Reads/MT.R{1,2}.fq.gz` (absolute paths are safest).

```bash
snakemake --snakefile Snakefile --configfile config/config_test.yaml \
          --use-conda --conda-prefix ./conda --cores 1 --dry-run
snakemake --snakefile Snakefile --configfile config/config_test.yaml \
          --use-conda --conda-prefix ./conda --cores 28 --printshellcmds
```

The first run also builds the remaining conda environments. On a Slurm cluster use the profile instead of
`--cores`: see [Running IMP3](running.md#on-a-slurm-cluster).

## 5. Your own samples

Copy `config/config.imp.yaml`, set the [inputs](inputs.md) and the [configuration](configuration.md), and
[run](running.md) it.

## Third-party code (`submodules/`)

Code that is not on conda lives in `submodules/` of the IMP3 folder. It is not part of the repository (third-party
licences); it is fetched either by the install workflow (step 3) or by setup rules on the first run, which need
internet access on the node where they run.

| folder | content | installed by | needed for |
|---|---|---|---|
| `binny/` | binny (fork, branch `imp`), with its marker database | install workflow (`binny_setup`) | binner `binny` |
| `cdhit/` | CD-HIT, compiled for long sequences | install workflow (`cdhit_compile`) | assembly merging |
| `KrakenTools/` | `combine_kreports.py` | install workflow (`krakentools_download`) | Kraken2 reports |
| `comebin/`, `comebin.ready` | COMEBin fork (improved clustering); patches the COMEBin environment | first run (`setup_comebin_scripts_from_fork`) | binner `COMEBin` |
| `bakta/`, `bakta.ready` | Bakta fork (dynamic database fetch); patches the Bakta environment | first run (`setup_bakta_functions_from_fork`) | annotation |

The first run also leaves marker files in `workflow/envs/` (`qm.ready`, `motus.ready`, `genome_updater.ready`, ...):
setup steps that patched an environment or downloaded into it (quickmerge wrapper, mOTUs database). Delete a marker
to redo its step.

Check an installation:

```bash
ls submodules/   # binny  cdhit  KrakenTools  (after the first run also: bakta  comebin  *.ready)
```

A second IMP3 checkout of the same version (e.g. a new clone next to a working installation) can reuse all of it:

```bash
cp -a /path/to/old/IMP3/submodules /path/to/new/IMP3/
cp -n /path/to/old/IMP3/workflow/envs/*.ready /path/to/new/IMP3/workflow/envs/
```

The conda environments are reused as well if the new checkout uses the same `--conda-prefix` and its
`workflow/envs/*.yaml` are unchanged. With the [IMP3 container](portable.md#the-imp3-container), another site gets
`submodules/` as `imp3_submodules.tar.gz` next to the image.

## Snakemake 7 shadow patch

Replace `unshadow_output` in `snakemake/dag.py` of your Snakemake environment
(`$CONDA_PREFIX/lib/python3.*/site-packages/snakemake/dag.py`) with:

```python
    def unshadow_output(self, job, only_log=False):
        """Move files from shadow directory to real output paths."""
        if not job.shadow_dir or not job.expanded_output:
            return

        # Combine outputs and logs based on the 'only_log' flag
        if only_log:
            files = job.log
        else:
            files = list(chain(job.expanded_output, job.log))

        # Create a mapping from absolute paths to output objects
        output_path_to_obj = {os.path.abspath(output.file): output for output in files}
        output_paths = set(output_path_to_obj.keys())
        top_level_paths = set()

        # Identify top-level outputs (those not nested within other outputs)
        for path in output_paths:
            is_top_level = True
            for other_path in output_paths:
                if other_path == path:
                    continue
                if other_path.startswith(path + os.sep):
                    is_top_level = False
                    break
            if is_top_level:
                top_level_paths.add(path)

        logger.debug(f"Top-level outputs identified for unshadowing: {top_level_paths}")

        # Move each top-level output
        for real_output_path in top_level_paths:
            real_output = output_path_to_obj[real_output_path]
            shadow_output = job.shadowed_path(real_output).file

            # Remake absolute symlinks as relative
            if os.path.islink(shadow_output):
                dest = os.readlink(shadow_output)
                if os.path.isabs(dest):
                    rel_dest = os.path.relpath(dest, job.shadow_dir)
                    os.remove(shadow_output)
                    os.symlink(rel_dest, shadow_output)
                    logger.debug(f"Converted absolute symlink {shadow_output} to relative symlink {rel_dest}")

            # Skip if shadow and real paths are the same
            if os.path.realpath(shadow_output) == os.path.realpath(real_output.file):
                continue

            logger.debug(f"Moving shadow output {shadow_output} to destination {real_output.file}")

            try:
                # Ensure the parent directory exists
                os.makedirs(os.path.dirname(real_output.file), exist_ok=True)

                # Move the top-level output
                shutil.move(shadow_output, real_output.file)
                logger.debug(f"Successfully moved {shadow_output} to {real_output.file}")

            except FileNotFoundError:
                logger.error(f"Source file or directory not found: {shadow_output}")
                raise
            except Exception as e:
                logger.error(f"Error moving {shadow_output} to {real_output.file}: {e}")
                raise

        # Verify that all specified outputs have been moved
        missing_outputs = [output.file for output in files if not os.path.exists(output.file)]
        if missing_outputs:
            logger.error(f"The following output files were not successfully moved: {missing_outputs}")
            raise FileNotFoundError(f"Unshadowing failed for outputs: {missing_outputs}")
        else:
            logger.debug("All specified outputs have been successfully moved.")

        # Cleanup the shadow directory if it's empty
        try:
            shutil.rmtree(job.shadow_dir)
            logger.debug(f"Shadow directory {job.shadow_dir} has been removed.")
        except Exception as e:
            logger.error(f"Failed to remove shadow directory {job.shadow_dir}: {e}")
            raise
```
