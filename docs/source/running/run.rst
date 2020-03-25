.. _run_IMP:

============
Running IMP3
============

You can run IMP3 by calling Snakemake with the appropriate settings:

.. code-block:: console

  snakemake /path/to/IMP3/Snakefile --configfile config.file.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cores corenumber

All paths can be relative or absolute. ``--cores`` specifies the number of threads you want to maximally use at the same time in addition to the
thread running the main snakemake (but note :ref:`extra limitations <hardware>` set in the config file).

.. _run_report:

If you'd like to have a **report** which lists all conda environments, rules and running stats, in addition
to the visual outputs, run snakemake with the ``--report`` option:

.. code-block:: console

  snakemake /path/to/IMP3/Snakefile --configfile config.file.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cores corenumber --report

.. _run_cluster:

If you want to use a batch submission system to run IMP3, add the appropriate arguments (more details :ref:`here <cluster_config>):

.. code-block:: console

  snakemake /path/to/IMP3/Snakefile --configfile config.file.yaml --use-conda --conda-prefix /path/to/IMP3/conda --cores corenumber --cluster-config /path/to/IMP3/config/cluster.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu}{params.mem} {cluster.threads}{threads} {cluster.partition} {cluster.nodes}" 


