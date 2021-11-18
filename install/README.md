# About

A `snakemake` workflow to setup all requirements for `IMP3`.

Config file: `config/config.install.yaml`

```bash
# from parent directory:
# activate the main conda environment
conda activate <your_imp3_conda_env>
# dry-run
snakemake -s install/Snakefile --configfile config/config.install.yaml --cores 1 -rpn
# execute (change the number of cores as needed)
snakemake -s install/Snakefile --configfile config/config.install.yaml --cores 1 -rp
```