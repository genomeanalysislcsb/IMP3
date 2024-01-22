# Prerequisites

## Conda

Conda is the recommended way of installing Snakemake, as it enables Snakemake to handle sofware dependencies of your workflow.

The recommended choice is Mambaforge which not only provides the required Python and Conda commands, but also includes Mamba an extremely fast and robust replacement for the Conda package manager which is highly recommended. The default conda solver is a bit slow and sometimes has issues with selecting the latest package releases. Therefore, we recommend to in any case use Mamba.

We use micromamba as it has all the packages necessary to run the IMP pipeline

```console
$ conda install -n base -c conda-forge micromamba
```

The databases used will take up a large amount of space so at least 1TB should be free for those databases. 

As the pipeline is quite resource-intensive, it is recommended to run it on an hpc if you have access to one. 
