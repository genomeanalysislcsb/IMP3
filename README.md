# IMP version 3

## Quickstart
Here is a quick guide on the installation and test run of IMP3. Please check out the longer description 
[here](https://imp3.readthedocs.io/en/latest/):

Note: All commands below are meant to be executed from within the imp3 directory after cloning. Cluster execution using 
profiles supports slurm only at the moment.

1) Clone the repository with git:
    ```shell
    git clone https://gitlab.lcsb.uni.lu/IMP/imp3.git -b binning_rework_v011
    cd imp3
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

    You need to choose a path to place the conda environments using Snakemakes `--conda-prefix` parameter, by default it
    will be put in `conda` in the IMP3 dir.

    Note that the installation might take a long time, if you selected large databases.
    ```shell
    # Activate environment
    conda activate my_snakemake_env
    
    # Set conda dir
    my_conda_dir="./conda"
    
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
              --cores 14 \
              --use-conda \
              --conda-prefix ${my_conda_dir} \
              --reason \
              --printshellcmds
    ```

5) Install remaining conda envs and run a test:
    
    You can use Snakemakes profile feature, to have it submit different jobs to the cluster or run it locally.
    
    First you will have to set up the test config file at `config/config.test.yaml`:
    1) Set up the specifications of your cluster nodes/ whatever kind of machine you have. If you don't have a separate
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
