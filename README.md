# IMP3

[![docs](https://img.shields.io/badge/docs-genomeanalysislcsb.github.io%2FIMP3-blue)](https://genomeanalysislcsb.github.io/IMP3/)

Integrated Meta-omic Pipeline, version 3: a reproducible Snakemake workflow for the integrated analysis of
metagenomic and metatranscriptomic data — preprocessing, hybrid assembly, annotation, read-based taxonomy, binning
into MAGs (CheckM2, GTDB-Tk) and per-sample statistics. Runs on metagenomes or metatranscriptomes alone, too.

**Documentation: <https://genomeanalysislcsb.github.io/IMP3/>**

- [Installation and test run](https://genomeanalysislcsb.github.io/IMP3/installation/)
- [Configuration](https://genomeanalysislcsb.github.io/IMP3/configuration/) ·
  [Running](https://genomeanalysislcsb.github.io/IMP3/running/) ·
  [Output files](https://genomeanalysislcsb.github.io/IMP3/outputs/)

```bash
git clone https://github.com/genomeanalysislcsb/IMP3.git && cd IMP3
conda env create -n snakemake_env -f requirements.yaml          # Snakemake 7.32
conda activate snakemake_env
snakemake -s install/Snakefile --configfile config/config.install.yaml --use-conda --cores 8   # tools + databases
snakemake --profile cluster_profile/ --configfile config/my_sample.yaml                        # one sample
```

> **Snakemake 7 note:** IMP3 needs a patched `unshadow_output` in Snakemake's `dag.py` for the shadow directive;
> see [Installation](https://genomeanalysislcsb.github.io/IMP3/installation/#snakemake-7-and-the-shadow-directive).

Third-party tools (binny, Bakta, COMEBin, CD-HIT, ESOM, KrakenTools) are installed by the install workflow and are
not part of this repository.

## Downstream

- [IMP explorer](https://github.com/genomeanalysislcsb/imp-explorer): study-level interactive exploration of IMP3 results
- [sRNA-IMP](https://github.com/genomeanalysislcsb/sRNA-IMP): small-RNA analysis on the IMP3 assembly

## History and licence

IMP3 continues IMP ([Narayanasamy et al. 2016](https://doi.org/10.1186/s13059-016-1116-8)). Development moved from
GitLab to GitHub with its full history. MIT licence (see `LICENSE`).
