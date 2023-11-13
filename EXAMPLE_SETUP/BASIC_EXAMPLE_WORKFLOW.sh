#!/bin/bash

##############################################################
##### DISCLAIMER #############################################
##############################################################
# Use this as a basic guide to run IMP3 on your data
# This script is not meant to be run as is, but rather to be used as a template
# You can copy this directory to your project directory and modify files accodring to your needs

##############################################################
##### Functions ##############################################
##############################################################
# Example usage
# run_parallel_snakemake true true 4 ${cluster_profile_aion} "/path/to/imp_dir" "/path/to/proj_dir" "log_prefix" "my_rule" "${aion_conf_array[@]}"
# To run full workflow, pass until_rule as empty string: ""
run_parallel_snakemake() {
    local parallel_dry_run=$1
    local snakemake_dry_run=$2
    local parallel_jobs=$3
    local cluster_profile=$4
    local imp_dir=$5
    local proj_dir=$6
    local log_prefix=$7
    local until_rule=$8
    local conf_array=("${@:9}")  # Get all arguments starting from the 8th one as an array
    
    # Log of parallel snakemake runs
    local main_log="${proj_dir}/run_${log_prefix}"

    # Determine if parallel dry-run
    local parallel_dry_run_cmd=""
    if [ "$parallel_dry_run" = "true" ]; then
        parallel_dry_run_cmd="--dry-run"
        main_log+="_parallel_dry_run"
    else
        main_log+=""
    fi

    # Determine if snakemake dry-run
    local snakemake_dry_run_cmd=""
    if [ "$snakemake_dry_run" = "true" ]; then
        snakemake_dry_run_cmd="--dry-run"
        main_log+="_snakemake_dry_run.log"
        log_prefix+="_dry_run"
    else
        main_log+=".log"
    fi

    # Include the --until option if until_rule is provided
    local until_cmd=""
    if [ -n "$until_rule" ]; then
        until_cmd="--until ${until_rule}"
    fi

    parallel -j ${parallel_jobs} ${parallel_dry_run_cmd} \
        "snakemake --profile ${cluster_profile} " \
        "--snakefile ${imp_dir}/Snakefile " \
        "--cluster-cancel scancel " \
        "--cluster-status ${imp_dir}/runscripts/status.py " \
        "--configfile {} --unlock; " \
        "snakemake --profile ${cluster_profile} " \
        "--snakefile ${imp_dir}/Snakefile " \
        "--cluster-cancel scancel " \
        "--cluster-status ${imp_dir}/runscripts/status.py " \
        "${until_cmd} " \
        "--verbose " \
        "--rerun-triggers mtime " \
        "${snakemake_dry_run_cmd} " \
        "--configfile {} > ${proj_dir}/logs/${log_prefix}_{/.}.log 2>&1" ::: "${conf_array[@]}" > "${main_log}" 2>&1
}

##############################################################
##### Define directories and initial configuration files #####
##############################################################
# IMP3 directory
imp_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/20230330_IMP3" ### USER_INPUT ###
# Cluster profile Aion
cluster_profile_aion="${imp_dir}/EXAMPLE_SETUP/cluster_profile_aion" ### USER_INPUT ###
# Cluster profile Iris
cluster_profile_iris="${imp_dir}/EXAMPLE_SETUP/cluster_profile_iris" ### USER_INPUT ###

# Project directory, for configs logs and imp outputs
proj_dir="/path/to/my/project/dir" ### USER_INPUT ###
# IMP3 sample output dir
imp_out_dir="${proj_dir}/imp3" ### USER_INPUT ###
# Path where to put the sample config files
conf_dir="${proj_dir}/imp_configs" ### USER_INPUT ###

# Path to the sample list file, tab-separated, with columns: MG1 MG2 MG_SE MT1 MT2 MT_SE
# If a sample is missing a read type, represent it as `NA`
# !!! USE ABSOLUTE PATHS !!!
# Example: /mnt/isilon/projects/ecosystem_biology/local_tools/20230330_IMP3/EXAMPLE_SETUP/example_sample_read_file_list.tsv
# There is a helper script to generate this file: /mnt/isilon/projects/ecosystem_biology/local_tools/20230330_IMP3/EXAMPLE_SETUP/build_sample_read_path_table_helper.sh
# Adapt to your needs
sample_list="${proj_dir}/path/to/samples.tsv" ### USER_INPUT ###

# The dir this file, the examples and the helper scripts are in
script_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/20230330_IMP3/EXAMPLE_SETUP" ### USER_INPUT ###

##############################################################
##### Make sample configs ####################################
##############################################################

# !!! CRITICAL !!! Before you create configs from the templates, make sure you have the correct parameters for your project set in BOTH template files !!! CRITICAL !!!
# Templates: /mnt/isilon/projects/ecosystem_biology/local_tools/20230330_IMP3/EXAMPLE_SETUP/config.imp.<aion/iris>_template.yaml

# Usage: 
# ./make_configs.sh [IMP_DIR] [CONF_TEMPLATE] [CONF_DIR] [SAMPLE_LIST] [HAS_HEADER=true] [IMP_RUN_OUT_DIR=./]
#
# Expects empty fields to be `NA`
#
# Parameters:
# IMP_DIR: Import directory path.
# CONF_TEMPLATE: Configuration template path.
# CONF_DIR: Output config directory path.
# SAMPLE_LIST: Sample list file path.
# HAS_HEADER: Sample list has header (default: true).
# IMP_RUN_OUT_DIR: Output directory for imp runs (default: current directory).

# Make out dir
mkdir -p ${imp_out_dir}

# !!! ENSURE CORRECT PARAMS FOR PROJECT !!!
for cluster in "aion" "iris"; do
  echo "Creating configs for cluster: ${cluster}"
  mkdir -p ${conf_dir}/${cluster}
  ${script_dir}/make_configs.sh \
      ${imp_dir} \
      ${script_dir}/config.imp.${cluster}_template.yaml \
      ${conf_dir}/${cluster} \
      ${sample_list}\
      true \
      ${imp_out_dir}
done

##############################################################
##### Make sample configs ####################################
##############################################################

# !!! CRITICAL !!! Before you create configs from the templates, make sure you have the correct parameters for your project set in BOTH template files !!! CRITICAL !!!
# Templates: /mnt/isilon/projects/ecosystem_biology/local_tools/20230330_IMP3/EXAMPLE_SETUP/config.imp.<aion/iris>_template.yaml

# Usage: 
# ./make_configs.sh [IMP_DIR] [CONF_TEMPLATE] [CONF_DIR] [SAMPLE_LIST] [HAS_HEADER=true] [IMP_RUN_OUT_DIR=./]
#
# Expects empty fields to be `NA`
#
# Parameters:
# IMP_DIR: Import directory path.
# CONF_TEMPLATE: Configuration template path.
# CONF_DIR: Output config directory path.
# SAMPLE_LIST: Sample list file path.
# HAS_HEADER: Sample list has header (default: true).
# IMP_RUN_OUT_DIR: Output directory for imp runs (default: current directory).

# Make out dir
mkdir -p ${imp_out_dir}

# !!! ENSURE CORRECT PARAMS FOR PROJECT !!!
for cluster in "aion" "iris"; do
  echo "Creating configs for cluster: ${cluster}"
  mkdir -p ${conf_dir}/${cluster}
  ${script_dir}/make_configs.sh \
      ${imp_dir} \
      ${script_dir}/config.imp.${cluster}_template.yaml \
      ${conf_dir}/${cluster} \
      ${sample_list}\
      true \
      ${imp_out_dir}
done

##############################################################
##### Create config file arrays ##############################
##############################################################
# Ensure arrays are cleared
unset aion_conf_array
unset iris_conf_array
# Iterate over sample dirs in the imp output dir
# If the status file `all.done`` does not exist in the directory, 
# add its corresponding yaml configuration to the cluster conf_array
for cluster in "aion" "iris"; do
  declare -a "${cluster}_conf_array"
  for i in "${conf_dir}/${cluster}"/*; do 
    if [[ ! -f "${imp_out_dir}/${i}/status/all.done" ]]; then
      sample="${i##*/}"
      eval "${cluster}_conf_array+=(\"${conf_dir}/${cluster}/${sample}\")"
    fi
  done
done


echo "Samples to run IMP workflow on: ${aion_conf_array[@]}"

##############################################################
##### Activate Snakemake environment #########################
##############################################################
# Choose: conda mamba micromamba
env_manager="micromamba" ### USER_INPUT ###
s_env="/mnt/isilon/projects/ecosystem_biology/local_tools/envs/snakemake_env"
${env_manager} activate ${s_env}

# Dont forget to paste the `run_parallel_snakemake` function in your terminal prior to running the rest of the workflow

##############################################################
##### Make log dirs ##########################################
##############################################################
mkdir -p ${proj_dir}/logs

##############################################################
##### Define misc parameters #################################
##############################################################
# Number of IMP runs to run in parallel
parallel_jobs=10 ### USER_INPUT ###
# When param assemlby["merge"] is set to "quickmerge", use: quickmerge_join_h1h2_with_h2h1
# When param assemlby["merge"] is set to ""/None, use: skip_merge_assembly
stage_1_until_rule="quickmerge_join_h1h2_with_h2h1" ### USER_INPUT ###

##############################################################
##### Run stage 1 on Aion ####################################
##############################################################
# Run steps: Preprocessing, Assembly, etc up to before Kraken tax annotations

# Dry-run  (snakemake, not parallel)
# The main log will be in ${proj_dir}/logs/run_stage_1_snakemake_dry_run.log
# Check the sample logs at ${proj_dir}/logs/run_stage_1_snakemake_dry_run_<sample_name>.log for errors
run_parallel_snakemake false true ${parallel_jobs} ${cluster_profile_aion} ${imp_dir} ${proj_dir} "stage_1" "quickmerge_join_h1h2_with_h2h1" "${aion_conf_array[@]}"
# Run stage 1
# This might take a long time, depending on the number of samples and the cluster load
# Ensure the command is running in a persistent screen/tmux session
# The log will be in ${proj_dir}/logs/run_stage_1.log
run_parallel_snakemake false false ${parallel_jobs} ${cluster_profile_aion} ${imp_dir} ${proj_dir} "stage_1" "quickmerge_join_h1h2_with_h2h1" "${aion_conf_array[@]}"

##############################################################
##### Ensure all samples have assembly files #################
##############################################################
# Check all sample directories for assembly files, except for tmp dir

for i in ${imp_out_dir}/*; do
  if [[ "${i}" != "${imp_out_dir}/tmp" ]]; then
    if ! ls ${i}/Assembly/*.assembly.merged.fa 1> /dev/null 2>&1; then
      echo "WARNING: No assembly files found for ${i}. Check IMP logs for errors."
    fi
  fi
done

##############################################################
##### Run stage 2 on Iris ####################################
##############################################################
# Run steps: Kraken tax annotations on reads and contigs

seq_type='reads'

# Dry-run (snakemake, not parallel) reads
run_parallel_snakemake false true ${parallel_jobs} ${cluster_profile_iris} ${imp_dir} ${proj_dir} "stage_2_${seq_type}" "kraken_${seq_type}" "${iris_conf_array[@]}"
# Run stage 2 reads
run_parallel_snakemake false false ${parallel_jobs} ${cluster_profile_iris} ${imp_dir} ${proj_dir} "stage_2_${seq_type}" "kraken_${seq_type}" "${iris_conf_array[@]}"

seq_type='contigs'
# Dry-run (snakemake, not parallel) contigs
run_parallel_snakemake false true ${parallel_jobs} ${cluster_profile_iris} ${imp_dir} ${proj_dir} "stage_2_${seq_type}" "kraken_${seq_type}" "${iris_conf_array[@]}"
# Run stage 2 contigs
run_parallel_snakemake false false ${parallel_jobs} ${cluster_profile_iris} ${imp_dir} ${proj_dir} "stage_2_${seq_type}" "kraken_${seq_type}" "${iris_conf_array[@]}"

##############################################################
##### Ensure all samples have kraken read and contig files ###
##############################################################
for i in ${imp_out_dir}/*; do
  if [[ "${i}" != "${imp_out_dir}/tmp" ]]; then
    if ! ls ${i}/Analysis/taxonomy/kraken/m[g,t].kraken.report 1> /dev/null 2>&1; then
        echo "WARNING: No Kraken read annotation files found for ${i}. Check IMP logs for errors."
    fi
    if ! ls ${i}/Analysis/taxonomy/kraken/m*.kraken.contig.report 1> /dev/null 2>&1; then
        echo "WARNING: No Kraken contig annotation files found for ${i}. Check IMP logs for errors."
    fi
  fi
done

##############################################################
##### Run stage 3 on Aion ####################################
##############################################################
# Run steps: Binning, Annotations, Stats, etc

# Dry-run (snakemake, not parallel)
run_parallel_snakemake false true ${parallel_jobs} ${cluster_profile_aion} ${imp_dir} ${proj_dir} "stage_3" "" "${aion_conf_array[@]}"
# Run stage 3
run_parallel_snakemake false false ${parallel_jobs} ${cluster_profile_aion} ${imp_dir} ${proj_dir} "stage_3" "" "${aion_conf_array[@]}"

##############################################################
##### Ensure all samples runs are complete ###################
##############################################################
for i in ${imp_out_dir}/*; do
  if [[ "${i}" != "${imp_out_dir}/tmp" ]]; then
    if ! ls ${i}/status/all.done 1> /dev/null 2>&1; then
        echo "WARNING: IMP run for ${i} is not complete. Check IMP logs for errors."
    fi
  fi
done

##############################################################
##### Additional commands ####################################
##############################################################

# Just unlock
sample="sample_name"
snakemake --profile ${cluster_profile_aion} \
          --snakefile ${imp_dir}/Snakefile  \
          --unlock \
          --verbose  \
          --configfile ${conf_dir}/${sample}.yaml

# Single run
sample="sample_name"
until_rule="quickmerge_join_h1h2_with_h2h1"
snakemake --profile ${cluster_profile_aion}  \
          --snakefile ${imp_dir}/Snakefile  \
          --cluster-cancel scancel  \
          --cluster-status ${imp_dir}/runscripts/status.py \
          --until ${until_rule} \
          --verbose  \
          --configfile ${conf_dir}/${sample}.yaml > ${proj_dir}/logs/${sample}_single_run_01.log 2>&1

# Basic local run. No profile
sample="sample_name"
conda_env="/mnt/isilon/projects/ecosystem_biology/local_tools/envs/imp3"
max_core
snakemake --snakefile ${imp_dir}/Snakefile  \
          --keep-going True \
          --rerun-incomplete True \
          --reason True \
          --printshellcmds True \
          --use-conda True \
          --conda-prefix ${conda_env} \
          --verbose  \
          --configfile ${conf_dir}/${sample}.yaml > ${proj_dir}/logs/${sample}_single_run_local_01.log 2>&1
