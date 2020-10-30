#!/bin/bash -l

# slurm settings if executed by sbatch; adjust to your needs
#SBATCH -J IMP3_Launcher
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=7-00:00:00
#SBATCH -p batch
#SBATCH -q long

# ARGS: positional arguments
#   1) snakemake env. (name OR path)
#   2) snakemake config file (see config/config.imp.yaml)
#   3) snakemake slurm config file (see config/slurm.config.yaml)

SMK_ENV=$1 # conda activate will fail if not specified correctly
SMK_CONFIG=$2
SMK_SLURM=$3
SMK_JOBS=5
SMK_CLUSTER="sbatch --partition {cluster.partition} --qos {cluster.qos} {cluster.nodes} \
--mem-per-cpu {params.mem} -n {threads} -t {params.runtime} --job-name={cluster.job-name}"

# Checks
echo "conda activate: ${SMK_ENV}"
conda activate ${SMK_ENV}

for ff in ${SMK_CONFIG} ${SMK_SLURM}; do
    if [[ -f "${ff}" && -s "${ff}" ]]; then 
        echo "Found ${ff}"
    else 
        echo "Error: ${ff} does not exist or is empty"
        exit 1
    fi
done

# Pipeline
snakemake -s Snakefile -rp --jobs ${SMK_JOBS} --local-cores 1 \
--configfile ${SMK_CONFIG} --use-conda \
--cluster-config ${SMK_SLURM} --cluster "${SMK_CLUSTER}"