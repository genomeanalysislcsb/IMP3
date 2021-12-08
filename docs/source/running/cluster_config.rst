.. _cluster_config:

=====================
Cluster configuration
=====================

`Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ workflows, and hence IMP3, can be run in
`cluster mode <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution>`_, i.e. `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ 
will submit all the rules (except if explicitly stated) to a batch submission system. For `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ 
to know how to submit to the cluster, the command and arguments for submission can be stated in the ``snakemake`` call (here for a `slurm <https://slurm.schedmd.com/documentation.html>`_ 
batch scheduling system):

.. code-block:: console

  snakemake -s /path/to/IMP3/Snakefile --configfile /path/to/sample.config.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cluster-config /path/to/IMP3/config/slurm.config.yaml --cluster "sbatch -t{params.runtime} --mem-per-cpu {params.mem} -n {threads} -p batch" 

The parameters (given in **brackets {}**) are predefined and will be taken from the ``params`` section of each `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ rule, 
so the user do not need to set them for every workflow task. More generally, a **cluster configuration file** (e.g. ``slurm.config.yaml``) can be used. Currently, 
IMP3 currently supplies cluster config files for `slurm <https://slurm.schedmd.com/documentation.html>`_ , for `PBS <https://www.openpbs.org/>`_ and for the `(Oracle) grid engine <https://en.wikipedia.org/wiki/Oracle_Grid_Engine>`_. 

If other scheduling systems are used, the **cluster configuration file** has to be adjusted. For help, please contact :ref:`IMP support <imp_support>`.

.. code-block:: console

  snakemake -s /path/to/IMP3/Snakefile --configfile /path/to/sample.config.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cluster-config /path/to/IMP3/config/slurm.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu}{params.mem} {cluster.threads}{threads} {cluster.partition} {cluster.nodes}" 

A short description for the IMP3 **cluster configuration file** is given below:

.. code-block:: yaml

  __default__:
	call: "sbatch"
	nodes: ""
	mem_per_cpu: "--mem-per-cpu "
	partition: "-p batch"
	quality: "-qos qos-batch"
	runtime: "-t"
	threads: "-n"
  
  mg_filtering:
	nodes: "-N 1"
	partition: "-p bigmem"
	quality: "-qos qos-bigmem"
	
	
The **cluster configuration file** is in `yaml <https://en.wikipedia.org/wiki/YAML>`_ format. There is a section about the default submit command at the top.
In addition, the steps that require more RAM (``bigmem``) are listed, together with the arguments to achieve ``bigmem`` submissions.
The user needs to fill in the command for submission and the way arguments are stated - we suggest using one of the
provided cluster config files (``config/slurm.config.yaml``, ``config/slurm_simple.config.yaml``, ``config/pbs.config.yaml``, or ``config/sge.config.yaml``) as template.
Please consider submitting a pull request on the `IMP3 gitlab repository <https://git-r3lab.uni.lu/IMP/IMP/>`_ for new cluster config files.
