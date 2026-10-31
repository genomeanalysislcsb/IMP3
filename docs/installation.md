# Installation and test run
A guide to installing IMP3 with its databases and running the test sample.

!!! note "Third-party tools are installed, not bundled"
    binny, Bakta, COMEBin, CD-HIT, ESOM and KrakenTools are not part of this repository. The install workflow
    (`install/Snakefile`) clones or downloads them into `submodules/` and sets up the databases.

Note: All commands below are meant to be executed from within the imp3 directory after cloning. Cluster execution using 
profiles supports slurm only at the moment.

1) Clone the repository with git:
    ```shell
    git clone https://github.com/genomeanalysislcsb/IMP3.git
    cd IMP3
    ```

2) Create a Snakemake environment or add some necessary packages to your existing one
    ```shell
    # Create new one
    my_package_manager='conda'  # conda mamba micromamba
    ${my_package_manager} env create -n snakemake_env -f requirements.yaml
    ```
    or
    ```shell
    # Add to existing
    my_package_manager='conda'  # conda mamba micromamba
    ${my_package_manager} activate my_snakemake_env
    ${my_package_manager} install pigz=2.6 filetype=1.2.0 -c conda-forge
    ```

3) Set install configuration:

   You need to set some basic parameters in the config file `config/config.install.yaml` first:

   1) Set `db_path` to where you want to keep your databases. Make sure that there is enough space available as, 
      depending on which databases you choose for Kraken2 and Mantis, you will need up to ~1tb.

   2) Choose your desired Kraken2 database with the `kraken2` parameter, commenting out the unused ones

   3) Set the desired Mantis dbs, following the instructions in the config at `mantis`.
   
   4) Set `snakemake_env_path` under `binny` to the path of either your existing Snakemake environment, the one you
      might have created in step 2) (e.g. `/abs/path/to/imp3/conda/snakemake_env`) or set it to `in_path`, if it is
      always available to your system by default.

      (If you want to find the path to your Snakemake environment, you can use `echo ${CONDA_PREFIX}`,
      while the environment is active.)

4) Install databases and prep some tools:

    You need to choose a path to place the conda environments using Snakemakes `--conda-prefix` parameter, in the example it
    will be put in `conda` in the IMP3 dir.

    You might need to set your conda channel priority to flexible to be able to install some of the environments. To do
    so you can run `conda config --set channel_priority flexible` or do it manually in the file `~/.condarc`. It is
    also recommended to have either no or a basic default channel setup  in `~/.condarc` with only bioconda and conda-forge,
    e.g. like this:
    ```yaml
    channel_priority: flexible
    channels:
      - conda-forge
      - bioconda
      - anaconda
    ```

    Note that the installation might take a long time, if you selected large databases.
    ```shell
    # Activate environment (replace name with the one you chose)
    ${my_package_manager} activate snakemake_env
    
    # Set conda dir and amount of cores (at least 8)
    my_conda_dir="./conda"
    cores=14
    
    # Dry-run
    snakemake --snakefile install/Snakefile \
              --configfile config/config.install.yaml \
              --cores 1 \
              --use-conda \
              --conda-prefix ${my_conda_dir} \
              --reason \
              --printshellcmds \
              --verbose \
              --dry-run
              
    # Run install
    snakemake --snakefile install/Snakefile \
              --configfile config/config.install.yaml \
              --cores ${cores} \
              --use-conda \
              --conda-prefix ${my_conda_dir} \
              --reason \
              --printshellcmds
    ```

5) Install remaining conda envs and run a test:
    
    You can use Snakemakes profile feature, to have it submit different jobs to the cluster or run it locally.
    
    First you will have to set up the test config file at `config/config.test.yaml`:
    1) Set up the specifications of your cluster nodes/whatever kind of machine you have under `mem`. If you don't have a separate
       large memory system you can set all `big_mem` parameter to be identical to the normal ones.

    2) Set up all paths for `tmp_dir`, `outputdir`, `krakendb`, the paths at `mantis`, and `Metagenomics` 
    (test files are in the IMP3 dir at `test/Reads`), matching what was done during the installation.

    Optional:
    
    3) Set up `cluster_profile/config.yaml` paths for `configfile`, `default-resources`, `conda-prefix`, and
       `cluster-config`, as well as max numbers of cores to use and jobs to run in parallel.

    4) Dry-run:
        ```shell
        # Activate environment
        conda activate my_snakemake_env
        
        # Dry-run
        my_conda_dir="./conda"
        snakemake --snakefile Snakefile \
                  --configfile config/config.test.yaml \
                  --use-conda \
                  --conda-prefix ${my_conda_dir} \
                  --cores 1 \
                  --reason \
                  --printshellcmds \
                  --verbose \
                  --dry-run
        ```
   
    5) Regular (Ensure enough memory is available, e.g. for the large kraken db you will need ~450gb:
        ```shell
        # Set conda dir and run install
        my_conda_dir="./conda"
        # Set max cores
        n_cores=28
       
        snakemake --snakefile Snakefile \
                 --configfile config/config.test.yaml \
                 --cores ${n_cores} \
                 --use-conda \
                 --conda-prefix ${my_conda_dir} \
                 --reason \
                 --printshellcmds
        ```
        or profile:
        ```shell
        # Set conda dir and run install
        my_conda_dir="./conda"
        snakemake --profile cluster_profile/ \
                  --verbose
        ```

6) Set up your own project:
    1) Set up your config using `config/config.imp.yaml` as template the same way it was described before.
   
    Optional:
    
    2) Set up your profile using `cluster_profile/config.yaml` as template the same way it was described before.

## Snakemake 7 and the shadow directive

For IMP3 to work with the shadow directive (Snakemake v7.*) we modify `unshadow_output` in `dag.py` so it properly copies back files:
```
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
