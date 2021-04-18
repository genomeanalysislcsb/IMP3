.. _install:

===================
How to install IMP3
===================

To install IMP3, you need `conda <https://www.anaconda.com/>`_.

1. Clone this repository to your disk:

.. code-block:: console

  git clone --single-branch --branch master https://git-r3lab.uni.lu/IMP/imp3 ./IMP3

Change into the IMP3 directory:

.. code-block:: console

  cd IMP3


At this point, you have all the scripts you need to run the workflow using snakemake. You still need to download a number of databases, though (see point 4).

2. **maybe optional** Install snakemake via conda:
If you don't have a recent version of snakemake in your path, you can install snakemake in a conda environment, as `recommended by Snakemake <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_ :

.. code-block:: console

  conda install -c conda-forge mamba
  mamba create --prefix $PWD/conda/snakemake_env
  conda activate $PWD/conda/snakemake_env
  mamba install -c conda-forge -c bioconda snakemake
  conda deactivate


3. **optional, but highly recommended**: Set permissions / PATH:
IMP3 is meant to be used by multiple users. Set the permissions accordingly. I'd suggest:

- to have read access for all files for the users, **plus**:
- read, write and execution rights for the conda subfolder
- It can also be useful to make the default config.imp.yaml not writable once you've set the paths to the databases you want to use (see below).


4. Download databases:
*note*: if you already have some of these databases, you don't need to download them again. You need to link them into the directory that keeps all of the databases, though. The path to this directory has to be set in the file config/config.imp.yaml.

.. code-block:: yaml

  db_path: "/absolute/path/to/your/databases" 
 
You may also not need all of these databases, if you don't run the steps that use them. Here, we list recommended databases for every step:
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


5. Initialize conda environments:
This sets up the conda environments that will be usable by all users and will download more databases:

.. code-block:: console

  snakemake --cores 1 -s Snakefile --conda-create-envs-only --use-conda --conda-prefix ./conda --configfile config/config.imp_init.yaml --local-cores 1

This step will take several minutes to an hour.
I strongly suggest to **remove one line from the activation script** after the installation, namely the one reading: `R CMD javareconf > /dev/null 2>&1 || true`, because you don't need this line later and if two users run this at the same time it can cause trouble. You can do this by running:

.. code-block:: console

  sed -i "s/R CMD javareconf/#R CMD javareconf/" conda/*/etc/conda/activate.d/activate-r-base.sh



