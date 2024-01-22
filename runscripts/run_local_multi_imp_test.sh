# Run IMP3 test on local machine

# Set up variables
se_isilon="/mnt/isilon/projects/ecosystem_biology"

imp_dir="${se_isilon}/local_tools/IMP3_devel"
imp_conf_file="${se_isilon}/local_tools/IMP3_devel/config/config_multi_imp_test.yaml"

max_cores=4

# Conda environments
sm_env="${se_isilon}/local_tools/envs/snakemake_env_dev"
imp_env="${se_isilon}/local_tools/envs/imp3"

# Chose conda environment manager
# `conda` or `mamba` or `micromamba`
env_manager="micromamba"

# Activate conda environment
micromamba activate ${sm_env}

# Run snakemake
snakemake --snakefile ${imp_dir}/Snakefile  \
          --reason \
          --printshellcmds \
          --use-conda \
          --conda-prefix ${imp_env} \
          --verbose  \
          --cores ${max_cores} \
          --configfile ${imp_conf_file}
