.. _install:

===================
How to install IMP3
===================

To install IMP3, you need `conda <https://www.anaconda.com/>`_.

1. Clone this repository to your disk:

.. code-block:: console

  git clone --recurse-submodules --single-branch --branch master https://git-r3lab.uni.lu/IMP/imp3.git ./IMP3

Change into the IMP3 directory:

.. code-block:: console

  cd IMP3


At this point, you have all the scripts you need to run the workflow.
You still need to download a number of databases, though (see point 3).

If you want to use the **comfortable IMP3 wrapper**, follow the points 4-8.
If you don't want to use the wrapper, point 9 has a hint on the initialization of the conda environments. 

2. Create the main ``conda`` environment:

If you don't have ``conda``, follow `these instructions first <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.

Create the main environment to be able to execute ``snakemake`` workflows.

.. code-block:: console

  conda env create -f requirements.yaml --prefix ./conda/snakemake_env

To run any ``snakemake`` workflow, you have to activate the created environment:

.. code-block:: console

  conda activate ./conda/snakemake_env

3. Download databases:

*Note*: if you already have some of these databases, you don't need to download them again.
You need to link them into the directory that keeps all of the databases, though.
You may also not need all of these databases, if you don't run the steps that use them.

There is an installation workflow (see directory ``install/``) which will setup the requirements defined in a config file.
A default configuration is provided in the file `config/config.install.yaml`.

Here is an example how to execute the setup:

.. code-block:: console

  # activate the snakemake environment
  conda activate ./conda/snakemake_env
  # dry-run
  snakemake -s install/Snakefile --configfile config/config.install.yaml --cores 1 -rpn
  # execute (use conda, more cores)
  snakemake -s install/Snakefile --configfile config/config.install.yaml --cores 8 --use-conda --conda-prefix ./conda -rp
  conda deactivate

Please note that the cores, memory and runtime requirements depend on the used config file.

For reference, see below for more information on how the requirements have to be set up for the different workflow steps.

**Step preprocessing**:

- Trimmomatic adapters: adapters should saved in a subdirectory in your database directory called ``adapters/``.
  You can find the latest version of Trimmomatic adapters `here <https://github.com/usadellab/Trimmomatic/tree/main/adapters/>`_ .
  .. We also `supply <https://webdav-r3lab.uni.lu/public/R3lab/IMP/Trimmomatic-Src-0.32.zip>`_ such a directory (from Trimmomatic v.0.32) that you'd just need to unzip.
  You may also want to use different adapters, which you can just supply in FASTA format but make sure to use the format expected by Trimmomatic.
  .. You specify the adapter set in the config file.
- Genomes to filter your reads against: reference FASTA files (with extention ``*.fa``) should be saved in a subdirectory in your database directory called ``filtering/``.
  Usually, you want to remove phiX174 and host genome sequences.
- SortMeRNA databases: when processing metatranscriptomic reads, the rRNA reads have to be removed.
  These are filtered by SortMeRNA against the provided reference FASTA files.
  .. We supply these `here <https://webdav-r3lab.uni.lu/public/R3lab/IMP/sortmerna.2.0.tgz>`_ .
  .. You need to unzip this into your database directory, or manually download the sequences from https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases and put them into a subdirectory in your database directory called ``sortmerna``.

.. In the config file, you can specify multiple filtering files like so:

.. .. code-block:: yaml

..   filtering:
..     filter: "phiX174 hg38"

**Step analysis**:

*TODO: Mantis*

.. - hmms: IMP3 assigns gene functions using HMMer3. For this, you need a subdirectory called ``hmm`` in your database directory. Each collection of HMMs you want to include must sit in its own subdirectory of hmm, which you specify in the config file, e.g. hmm_DBs: "Resfams essential". You can use multiple collections of HMMs, of course. Here are some suggestions:

.. - `Resfams <http://www.dantaslab.org/resfams/>`_
.. - `essential genes <https://webdav-r3lab.uni.lu/public/R3lab/IMP/essential.hmm>`_
.. - `Pfam-A <ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz>`_ *note*: please name the folder holding these ``Pfam`` or ``Pfam_A`` (avoid minus in the filename)
.. - `KEGG <ftp://ftp.genome.jp/pub/db/kofam/>`_
.. Contact us, if you're interested in more in-house databases.

**Step taxonomy**:

- Kraken2: to run taxonomic classification with Kraken2, a directory containing the required databases files has to be put (or linked) in your database directory.
  .. You will need a kraken2 database in your database directory (or a link to one).
  See `here <https://benlangmead.github.io/aws-indexes/k2>`_ for a collection of Kraken2 databases.
  You can of course also use a custom database.
  .. You may also want to build your own. You just supply the name of the k2 database in the config file later, e.g.

.. .. code-block:: yaml

..   krakendb: "minikraken2"

- GTDB-tk: if you have performed **binning**, you need the `GTDB database <https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/>`_ (file ``gtdbtk_data.tar.gz``).
  .. Untar it it to a subdirectory called ``GTDB_tk`` within your database directory.

4. Adjust the file VARIABLE_CONFIG to your requirements (have a tab between the variable name and your setting):

- SNAKEMAKE_VIA_CONDA - set this to true, if you don't have snakemake in your path and want to install it via conda (recommended, because that way you have a current version). Leave empty, if you don't need an additional snakemake.
- LOADING_MODULES - insert a bash command to load modules, if you need them to run conda. Leave empty, if you don't need to load a module.
- SUBMIT_COMMAND - insert the bash command you'll usually use to submit a job to your cluster to run on a single cpu for the entirety of the run (a few hours to days). You can also include extra arguments, in particularly, if you need to bind jobs to one node (e.g. sbatch --nodelist=). You only need this whole setting, if you want to have the snakemake top instance running in a submitted job. You alternatively have the option to run it on the frontend via tmux. Leave empty, if you want to use this version and have `tmux <https://github.com/tmux/tmux/wiki/>`_ installed.
- SCHEDULER - insert the name of the scheduler you want to use (currently `slurm` or `sge`). This determines the cluster config given to snakemake, e.g. the cluster config file for slurm is config/slurm.config.yaml or config/slurm_simple.config.yaml. Also check that the settings in this file is correct. If you have a different system, `contact us <https://git-r3lab.uni.lu/IMP/imp3/-/issues/>`_ and feel free to submit new scheduler files.
- MAX_THREADS - set this to the maximum number of cores you want to be using in a run. If you don't set this, the default will be 50. Users can override this setting at runtime.
- BIGMEM_CORES - set this to the maximum number of bigmem cores you want to be using in a run. If you don't set this, the default will be 5.
- NORMAL_MEM_EACH - set the size of the RAM of one core of your normal copute nodes (e.g. 8G). If you're not planning to use binny to submit to a cluster, you don't need to set this.
- BIGMEM_MEM_EACH - set the size of the RAM of one core of your bigmem (or highmem) compute nodes. If you're not planning to use binny to submit to a cluster or don't have separate bigmem nodes, you don't need to set this.
- NODENAME_VAR - if you submit your workflow to a cluster that requires binding to a fixed node, set this variable according to the environmental storing the nodename (e.g. SLURMD_NODENAME for most slurm systems)
- BIND_JOBS_TO_MAIN - if you submit your workflow to a cluster that requires binding to a fixed node, set this variable to true

5. Decide how you want to run IMP3, if you let it submit jobs to the cluster:

Only do one of the two:

- if you want to submit the process running snakemake to the cluster:

.. code-block:: console

  cp runscripts/runIMP3_submit.sh runIMP3
  chmod 755 runIMP3

- if you want to keep the process running snakemake on the frontend using tmux:

.. code-block:: console

  cp runscripts/runIMP3_tmux.sh runIMP3
  chmod 755 runIMP3

6. **optional, but highly recommended**: Install snakemake (and other dependencies) via conda:

If you want to use snakemake via conda (and you've set SNAKEMAKE_VIA_CONDA to true), install the environment, as `recommended by Snakemake <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_ :
Please note that the IMP3 pipeline is tested only using package versions specified in `requirements.yaml`. Using different versions than the specified ones might cause issues.

.. code-block:: console

  conda env create -f requirements.yaml --prefix $PWD/conda/snakemake_env


7. **optional, but highly recommended**: Set permissions / PATH:

IMP3 is meant to be used by multiple users. Set the permissions accordingly. I'd suggest:

- to have read access for all files for the users, **plus**:
- execution rights for the runIMP3 file and the .sh scripts in the subfolder submit_scripts
- read, write and execution rights for the conda subfolder
- to add the IMP3 directory to your path.
- It can also be useful to make the VARIABLE_CONFIG file not-writable, because you will always need it. The same goes for config.imp.yaml once you've set the paths to the databases you want to use (see below).

8. Initialize conda environments:

This run sets up the conda environments that will be usable by all users and will download a database:

.. code-block:: console

  ./runIMP3 -i config/config.imp_init.yaml


This step will take several minutes to an hour. 

9. Initialize the conda environments without wrapper:

This sets up the conda environments that will be usable by all users and will download more databases:

.. code-block:: console

  snakemake --cores 1 -s Snakefile --conda-create-envs-only --use-conda --conda-prefix `pwd`/conda --configfile config/config.imp_init.yaml --local-cores 1

This step will take several minutes to an hour.
I strongly suggest to **remove one line from the activation script** after the installation, namely the one reading: `R CMD javareconf > /dev/null 2>&1 || true`, because you don't need this line later and if two users run this at the same time it can cause trouble. You can do this by running:

.. code-block:: console

  sed -i "s/R CMD javareconf/#R CMD javareconf/" conda/*/etc/conda/activate.d/activate-r-base.sh
