.. _install:

===================
How to install IMP3
===================

To install IMP3, you need `conda <https://www.anaconda.com/>`_.

1. Clone this repository to your disk:

.. code-block:: console

  git clone --single-branch --branch master https://git-r3lab.uni.lu/IMP/imp3.git ./IMP3

Change into the IMP3 directory:

.. code-block:: console

  cd IMP3


At this point, you have all the scripts you need to run the workflow using snakemake. You still need to download a number of databases, though (see point 2). If you want to use the **comfortable IMP3 wrapper**, follow the points 3-7. If you don't want to use the wrapper, point 8 has a hint on the initialization of the conda environments. 

2. Download databases:
*note*: if you already have some of these databases, you don't need to download them again. You need to link them into the directory that keeps all of the databases, though. You may also not need all of these databases, if you don't run the steps that use them. Here, we list recommended databases for every step:
**preprocessing**:

- Trimmomatic adapters: You need a subdirectory in your database directory called adapters that hold the adapters you'd like to trim. You can find the latest versions `here <https://github.com/usadellab/Trimmomatic/tree/main/adapters/>`_ . We also `supply <https://webdav-r3lab.uni.lu/public/R3lab/IMP/Trimmomatic-Src-0.32.zip>`_ such a directory (from Trimmomatic v.0.32) that you'd just need to unzip. You may also want to use different adapters, which you can just supply in fasta format. You specify the adapter set in the config file.
- genomes to filter your reads against: You should remove phiX174 sequences, as they are more than likely to come from the sequencing. You need to download this file and place it in a subdirectory of your database directory. In the config file, you specify the filtering file like so:

.. code-block:: yaml

  filtering:
    filter: "phiX174"

You may in addition want to remove reads mapping to a host genome. Just like the phiX174 genome, you add the genome file in fasta format (with file ending .fa) in the filtering directory. You can specify multiple filtering files like so:

.. code-block:: yaml

  filtering:
    filter: "phiX174 hg38"

- sortmeRNA databases: If you run IMP3 on metatranscriptomic reads, the rRNA reads are filtered by sortmeRNA. We supply these `here <https://webdav-r3lab.uni.lu/public/R3lab/IMP/sortmerna.2.0.tgz>`_ . You need to unzip this into your database directory, or manually download the sequences from https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases and put them into a subdirectory in your database directory called sortmerna.

**analysis**:

- hmms: IMP3 assigns gene functions using HMMer3. For this, you need a subdirectory called hmm in your database directory. Each collection of HMMs you want to include must sit in its own subdirectory of hmm, which you specify in the config file, e.g. hmm_DBs: "Resfams essential". You can use multiple collections of HMMs, of course. Here are some suggestions:

- `Resfams <http://www.dantaslab.org/resfams/>`_
- `essential genes <https://webdav-r3lab.uni.lu/public/R3lab/IMP/essential.hmm>`_
- `Pfam-A <ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz>`_ *note*: please name the folder holding these Pfam or Pfam_A (avoid minus in the filename)
- `KEGG <ftp://ftp.genome.jp/pub/db/kofam/>`_
Contact us, if you're interested in more in-house databases.

**taxonomy**:

- kraken: You will need a kraken2 database in your database directory (or a link to one). You can get many of those `here <https://benlangmead.github.io/aws-indexes/k2/>`_ . You may also want to build your own. You just supply the name of the k2 database in the config file later, e.g.

.. code-block:: yaml

  krakendb: minikraken2"

- GTDB_tk: If you have performed **binning**, you need the `GTDB database <https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/>`_ . Store it in a subdirectory called GTDB_tk within your database directory.


3. Adjust the file VARIABLE_CONFIG to your requirements (have a tab between the variable name and your setting):

- SNAKEMAKE_VIA_CONDA - set this to true, if you don't have snakemake in your path and want to install it via conda (recommended, because that way you have a current version). Leave empty, if you don't need an additional snakemake.
- LOADING_MODULES - insert a bash command to load modules, if you need them to run conda. Leave empty, if you don't need to load a module.
- SUBMIT_COMMAND - insert the bash command you'll usually use to submit a job to your cluster to run on a single cpu for a few days. You only need this, if you want to have the snakemake top instance running in a submitted job. You alternatively have the option to run it on the frontend via tmux. Leave empty, if you want to use this version and have `tmux <https://github.com/tmux/tmux/wiki/>`_ installed.
- SCHEDULER - insert the name of the scheduler you want to use (currently `slurm` or `sge`). This determines the cluster config given to snakemake, e.g. the cluster config file for slurm is config/slurm.config.yaml or config/slurm_simple.config.yaml. Also check that the settings in this file is correct. If you have a different system, `contact us <https://git-r3lab.uni.lu/IMP/imp3/-/issues/>`_ and feel free to submit new scheduler files.
- MAX_THREADS - set this to the maximum number of cores you want to be using in a run. If you don't set this, the default will be 50. Users can override this setting at runtime.
- BIGMEM_CORES - set this to the maximum number of bigmem cores you want to be using in a run. If you don't set this, the default will be 5.
- NORMAL_MEM_EACH - set the size of the RAM of one core of your normal copute nodes (e.g. 8G). If you're not planning to use binny to submit to a cluster, you don't need to set this.
- BIGMEM_MEM_EACH - set the size of the RAM of one core of your bigmem (or highmem) compute nodes. If you're not planning to use binny to submit to a cluster or don't have separate bigmem nodes, you don't need to set this.

4. Decide how you want to run IMP3, if you let it submit jobs to the cluster:
Only do one of the two:

- if you want to submit the process running snakemake to the cluster:

.. code-block:: console

  cp runscripts/runIMP3_submit.sh runIMP3
  chmod 755 runIMP3

- if you want to keep the process running snakemake on the frontend using tmux:

.. code-block:: console

  cp runscripts/runIMP3_tmux.sh runIMP3
  chmod 755 runIMP3


5. **optional, but recommended**: Install snakemake via conda:
If you want to use snakemake via conda (and you've set SNAKEMAKE_VIA_CONDA to true), install the environment, as `recommended by Snakemake <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_ :

.. code-block:: console

  conda install -c conda-forge mamba
  mamba create --prefix $PWD/conda/snakemake_env
  conda activate $PWD/conda/snakemake_env
  mamba install -c conda-forge -c bioconda snakemake
  conda deactivate


6. **optional, but highly recommended**: Set permissions / PATH:
IMP3 is meant to be used by multiple users. Set the permissions accordingly. I'd suggest:

- to have read access for all files for the users, **plus**:
- execution rights for the runIMP3 file and the .sh scripts in the subfolder submit_scripts
- read, write and execution rights for the conda subfolder
- to add the IMP3 directory to your path.
- It can also be useful to make the VARIABLE_CONFIG file not-writable, because you will always need it. The same goes for config.imp.yaml once you've set the paths to the databases you want to use (see below).

7. Initialize conda environments:
This run sets up the conda environments that will be usable by all users and will download a database:

.. code-block:: console

  ./runIMP3 -i config/config.imp_init.yaml


This step will take several minutes to an hour. 

8. Initialize the conda environments without wrapper:
This sets up the conda environments that will be usable by all users and will download more databases:

.. code-block:: console

  snakemake --cores 1 -s Snakefile --conda-create-envs-only --use-conda --conda-prefix `pwd`/conda --configfile config/config.imp_init.yaml --local-cores 1

This step will take several minutes to an hour.
I strongly suggest to **remove one line from the activation script** after the installation, namely the one reading: `R CMD javareconf > /dev/null 2>&1 || true`, because you don't need this line later and if two users run this at the same time it can cause trouble. You can do this by running:

.. code-block:: console

  sed -i "s/R CMD javareconf/#R CMD javareconf/" conda/*/etc/conda/activate.d/activate-r-base.sh
 


