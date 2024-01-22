###
sample="1226_1"

sample_dir="/mnt/isilon/projects/ecosystem_biology/MiBiPa/IMP/${sample}/run1"
test_dir="/mnt/lscratch/users/ohickl/expobiome/tests/srna_imp_test"
config_file="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel/config/config.imp.hybrid_test.yaml"

micromamba activate /mnt/isilon/projects/ecosystem_biology/local_tools/envs/snakemake_env_dev

imp_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"

snakemake --profile ${imp_dir}/cluster_profile_test_aion \
          --snakefile ${imp_dir}/Snakefile  \
          --cluster-cancel scancel  \
          --cluster-status ${imp_dir}/runscripts/status.py  \
          --verbose  \
          --configfile ${config_file} -n > ${imp_dir}/imp3_test_0001.log 2>&1
