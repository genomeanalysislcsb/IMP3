# About

This is a `snakemake` workflow to setup `IMP3` dependencies including required databases etc.

The requirements are defined in the config file: `config/config.install.yaml`.

Please note that the cores, memory and runtime requirements depend on the used config file.

```bash
# from parent directory:

# activate the main conda environment
conda activate <your_imp3_conda_env>

# dry-run
snakemake -s install/Snakefile --configfile config/config.install.yaml --cores 1 -rpn

# execute (use conda, more cores)
snakemake -s install/Snakefile --configfile config/config.install.yaml --cores 12 --use-conda -rp
```