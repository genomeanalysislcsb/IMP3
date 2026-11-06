# IMP3

[![docs](https://img.shields.io/badge/docs-genomeanalysislcsb.github.io%2FIMP3-blue)](https://genomeanalysislcsb.github.io/IMP3/)

Integrated Meta-omic Pipeline, version 3: a reproducible Snakemake workflow for the integrated analysis of
metagenomic and metatranscriptomic data — preprocessing, hybrid assembly, annotation, read-based taxonomy, binning
into MAGs (CheckM2, GTDB-Tk) and per-sample statistics. Runs on metagenomes or metatranscriptomes alone, too.

**Documentation: <https://genomeanalysislcsb.github.io/IMP3/>**

- [Installation and test run](https://genomeanalysislcsb.github.io/IMP3/installation/) ·
  [Databases](https://genomeanalysislcsb.github.io/IMP3/databases/)
- [Inputs](https://genomeanalysislcsb.github.io/IMP3/inputs/) ·
  [Configuration](https://genomeanalysislcsb.github.io/IMP3/configuration/) ·
  [Running](https://genomeanalysislcsb.github.io/IMP3/running/) ·
  [Many samples, any machine](https://genomeanalysislcsb.github.io/IMP3/portable/)
- [Workflow steps](https://genomeanalysislcsb.github.io/IMP3/steps/overview/) ·
  [Output files](https://genomeanalysislcsb.github.io/IMP3/outputs/) ·
  [FAQ](https://genomeanalysislcsb.github.io/IMP3/faq/)

```bash
git clone https://github.com/genomeanalysislcsb/IMP3.git && cd IMP3
conda env create -n snakemake_env -f requirements.yaml          # Snakemake 7.32
conda activate snakemake_env
snakemake -s install/Snakefile --configfile config/config.install.yaml --use-conda --conda-prefix conda --cores 8   # tools + databases
snakemake --profile cluster_profile/ --configfile config/my_sample.yaml                        # one sample
```

Many samples, on a Slurm cluster or a single machine, with the conda environments or the IMP3
Apptainer/Singularity container ([details](https://genomeanalysislcsb.github.io/IMP3/portable/)):

```bash
./imp3 init my_project --executor slurm --site generic --db /path/to/imp3_db   # or --executor local, --container imp3.sif
cd my_project                       # reads in samples/<sample>/ (file names with mg / mt)
/path/to/IMP3/imp3 run --dry-run && /path/to/IMP3/imp3 run   # --submit: masters as Slurm jobs instead of tmux
/path/to/IMP3/imp3 status
```

> **Snakemake 7 note:** IMP3 needs a patched `unshadow_output` in Snakemake's `dag.py` for the shadow directive;
> see [Installation](https://genomeanalysislcsb.github.io/IMP3/installation/#snakemake-7-and-the-shadow-directive).

Third-party tools (binny, Bakta, COMEBin, CD-HIT, KrakenTools) are installed by the install workflow and are
not part of this repository.

## Downstream

- [IMP explorer](https://github.com/genomeanalysislcsb/imp-explorer): study-level interactive exploration of IMP3 results
- [sRNA-IMP](https://github.com/genomeanalysislcsb/sRNA-IMP): small-RNA analysis on the IMP3 assembly

Known open points (removed or postponed features): [TODO.md](TODO.md).

The Apptainer image with all environments of a release is on Zenodo:
[doi:10.5281/zenodo.24003090](https://doi.org/10.5281/zenodo.24003090) (see [Many samples, any machine](https://genomeanalysislcsb.github.io/IMP3/portable/#download-the-image)).

When you publish results, cite IMP and the tools you ran: references in
[History and credits](https://genomeanalysislcsb.github.io/IMP3/about/#citing-the-tools).

## History and licence

IMP3 continues IMP ([Narayanasamy et al. 2016](https://doi.org/10.1186/s13059-016-1116-8)). Development moved from
GitLab to GitHub with its full history. MIT licence (see `LICENSE`).
