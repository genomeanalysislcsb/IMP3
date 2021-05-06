.. _run_IMP:

============
Running IMP3
============

You can run IMP3 by calling `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ with the appropriate settings:

.. code-block:: console

  snakemake -s /path/to/IMP3/Snakefile --configfile /path/to/config.file.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cores corenumber

``Snakefile`` defines your workflow that you want to execute.

``--configfile`` specifies the config file for your sample. More info you can find :ref:`here <configuration>`.

``--use-conda`` enables `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ to use `conda <https://docs.conda.io/en/latest/>`_ environments defined for your workflow.

``--conda-prefix`` path to your IMP3 conda directory for the IMP3 conda environments.

``--cores`` specifies the number of threads you want to maximally use at the same time in addition to the
thread running the main ``snakemake`` command (but note :ref:`extra limitations <hardware>` set in the config file).

All paths can be relative or absolute. 

.. _run_report:

If you'd like to have a **report** which lists all IMP3 conda environments, rules and running stats, in addition
to the visual outputs, run `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ with the ``--report`` option:

.. code-block:: console

  snakemake /path/to/IMP3/Snakefile --configfile config.file.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cores corenumber --report

.. _run_cluster:

If you want to use a batch submission system to run IMP3, add the appropriate arguments (more details :ref:`here <cluster_config>`):

.. code-block:: console

  snakemake /path/to/IMP3/Snakefile --configfile config.file.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cores corenumber --cluster-config /path/to/IMP3/config/cluster.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu}{params.mem} {cluster.threads}{threads} {cluster.partition} {cluster.nodes}" 


