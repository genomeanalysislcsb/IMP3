# IMP3

The Integrated Meta-Omics Pipeline, IMP, is a reproducible and modular pipeline for large-scale
integrated analysis of coupled metagenomic and metatranscriptomic data. IMP3 incorporates read
preprocessing, assembly (including optional iterative co-assembly of metagenomic and
metatranscriptomic data), analysis of microbial community structure (taxonomy) and function
(gene functions and abundances), binning, as well as summaries and visualizations.

IMP3 is implemented in [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html). It is user-friendly
and can be easily configured to perform single-omic metagenomic or metatranscriptomic analyses, as
well as single modules, e.g. only assembly or only binning.

IMP3 is being developed at the [Luxembourg Centre for Systems Biomedicine](https://wwwen.uni.lu/lcsb/research/eco_systems_biology) and the
[Swammerdam Institute for Life Sciences at the University of Amsterdam](https://www.uva.nl/en/profile/h/e/a.u.s.heintzbuschart/a.u.s.heintz-buschart.html).

```{toctree}
:caption: For Uni.lu users
:hidden: true
:maxdepth: 1
:name: unilu

unilu/unilu
```

```{toctree}
:caption: Installation
:hidden: true
:maxdepth: 1
:name: installation

installation/prerequisites
installation/quickstart
```

```{toctree}
:caption: IMP3
:hidden: true
:maxdepth: 1
:name: IMP3

installation/installation
running/run_idx
steps/Overview
output/Overview
```
```{toctree}
:caption: walkthrough
:hidden: true
:maxdepth: 1
:name: walkthrough

walkthrough/walkthrough
```

```{toctree}
:caption: Multi-IMP
:hidden: true
:maxdepth: 1
:name: MultiIMP

multiImp/multiimp
```

```{toctree}
:caption: FAQ
:hidden: true
:maxdepth: 1
:name: FAQ

FAQ
```



% * :ref:`genindex`

% * :ref:`modindex`

% * :ref:`search`
