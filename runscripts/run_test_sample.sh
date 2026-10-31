#!/bin/bash -l

# Extra args for snakemake
extra_args="-n" # --conda-create-envs-only

# Get current date and time
datetime=$(date '+%Y-%m-%d_%H-%M-%S')

# Run a test sample with IMP
env_manager="micromamba" # conda mamba micromamba

SM_ENV_DIR="/mnt/isilon/projects/ecosystem_biology/local_tools/envs/snakemake_env_dev"
${env_manager} activate ${SM_ENV_DIR}  #|| { echo "Error: Failed to activate env ${SM_ENV_DIR}"; exit 1; }

imp_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"
snakefile="${imp_dir}/Snakefile"

proj_path="/mnt/lscratch/users/ohickl/tests"
sample="imp_test_${datetime}"

log_path="${proj_path}/${sample}.log"

# Unlock to ensure no issues with snakemake
snakemake --snakefile ${snakefile} \
          --unlock \
          --quiet all

snakemake --profile ${imp_dir}/cluster_profile_test \
          --snakefile ${snakefile} \
          --verbose  \
          ${extra_args} >> ${log_path} 2>&1
