.. _cluster_config:

=====================
Cluster configuration
=====================

Snakemake workflows, and hence IMP3, can be run in
`cluster mode <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution>`_, i.e. Snakemake will submit all the
rules (except if explicitly stated) to a batch submission system. For Snakemake to know how to submit to the cluster,
the command and arguments for submission can be stated
in the `snakemake` call:

.. code-block:: console

  snakemake -s /path/to/IMP3/Snakefile --configfile /path/to/sample.config.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cluster-config /path/to/IMP3/config/slurm.config.yaml --cluster "sbatch -t{params.runtime} --mem-per-cpu {params.mem} -n {threads} -p batch" 

The parameters will be taken from the params section of each rule, so you don't need to worry about them.
More generally, a cluster configuration file can be used. Currently, IMP3 supplies two cluster
config files, for slurm and for grid engines. If you use either system, you're in luck and you can just run IMP as is:

.. code-block:: console

  snakemake -s /path/to/IMP3/Snakefile --configfile /path/to/sample.config.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cluster-config /path/to/IMP3/config/slurm.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu}{params.mem} {cluster.threads}{threads} {cluster.partition} {cluster.nodes}" 

**If you want to write your own cluster config file, here's an explanation:**
Here's what the top of the config file for slurm looks like:

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
	
	
The cluster configuration file is in .yaml format.
There is a section about the default submit command at the top.
In addition, the steps that require more RAM (bigmem nodes) are listed, together with the arguments to achieve bigmem submissions.
All you need to do, is to fill in the command for submission and the way arguments are stated - we suggest using one of the
existing cluster config files (config/slurm.config.yaml or config/sge.config.yaml) as template.
Please consider submitting a pull request on `IMP3's gitlab <https://git-r3lab.uni.lu/IMP/IMP/>`_ repository for new cluster config files.
