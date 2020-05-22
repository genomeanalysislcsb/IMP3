#!/bin/bash -l

##############################
# SLURM
# NOTE: used for this script only, NOT for the snakemake call below

#SBATCH -J IMP_PF
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=1-00:00:00
#SBATCH -p batch
#SBATCH --qos=qos-batch

##############################
# SNAKEMAKE

# conda env name
IMP_ENV="IMP_PF"
# number of cores for snakemake
IMP_CORES=5
# IMP config file
IMP_CONFIG="config/config.imp.yaml" # USER INPUT REQUIRED
# slurm config file
IMP_SLURM="config/slurm.config.yaml"
# slurm cluster call
IMP_CLUSTER="{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu}{params.mem} {cluster.threads}{threads} {cluster.partition} {cluster.nodes} {cluster.quality}"

##############################
# IMP

# activate the env
conda activate ${IMP_ENV}

# run the pipeline
snakemake -s Snakefile -rp --cores ${IMP_CORES} --configfile ${IMP_CONFIG} \
--use-conda --conda-prefix ${CONDA_PREFIX}/IMP \
--cluster-config ${IMP_SLURM} --cluster "${IMP_CLUSTER}"
