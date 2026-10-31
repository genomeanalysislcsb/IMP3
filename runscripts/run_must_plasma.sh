

env_manager="micromamba" # conda mamba micromamba
micromamba activate /mnt/isilon/projects/ecosystem_biology/local_tools/envs/snakemake_env

# Define directories and initial configuration files
proj_dir="/mnt/lscratch/users/ohickl/MUST-totalRNA-plasma"
imp_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"
conf_dir="${proj_dir}/imp_configs"
conf_array=(${conf_dir}/*.yaml)


# Clear the configuration array
declare -a conf_array=()

# Iterate over directories in the sample dir *
# If the status file all.done does not exist in the directory, 
# add its corresponding yaml configuration to the conf_array
sample_dir="${proj_dir}/path/to/samples/*"
for i in ${proj_dir}/imp3/M*; do 
  if [[ ! -f ${i}/status/all.done ]]; then
    sample=${i##*/}
    conf_array+=(${conf_dir}/${sample}.yaml)
  fi
done


conf_array=(${conf_dir}/M4.3-V1B.yaml ${conf_dir}/M1.4-V2.yaml)

echo ${conf_array[@]}

parallel -j 50 --dry-run \
        "snakemake --profile ${imp_dir}/cluster_profile_test_aion " \
        "--snakefile ${imp_dir}/Snakefile " \
        "--cluster-cancel scancel " \
        "--cluster-status ${imp_dir}/runscripts/status.py " \
        "--configfile {} --unlock; " \
        "snakemake --profile ${imp_dir}/cluster_profile_test_aion " \
        "--snakefile ${imp_dir}/Snakefile " \
        "--cluster-cancel scancel " \
        "--cluster-status ${imp_dir}/runscripts/status.py " \
        "--verbose " \
        "--configfile {} > ${proj_dir}/imp3/logs/{/.}_08.log 2>&1" ::: "${conf_array[@]}" > run_must_plasma_08.log 2>&1




# unlock
sample="M4.3-V1B" # M1.4-V2 M4.3-V1B
snakemake --profile /mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel/cluster_profile_test_aion  \
          --snakefile /mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel/Snakefile  \
          --unlock \
          --verbose  \
          --configfile /mnt/lscratch/users/ohickl/MUST-totalRNA-plasma/imp_configs/${sample}.yaml

# single
snakemake --profile /mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel/cluster_profile_test_aion  \
          --snakefile /mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel/Snakefile  \
          --cluster-cancel scancel  \
          --cluster-status /mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel/runscripts/status.py  \
          --until gather_stats \
          --verbose  \
          --configfile /mnt/lscratch/users/ohickl/MUST-totalRNA-plasma/imp_configs/${sample}.yaml > /mnt/lscratch/users/ohickl/MUST-totalRNA-plasma/imp3/logs/${sample}_09.log 2>&1