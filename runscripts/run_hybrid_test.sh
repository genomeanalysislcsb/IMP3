
# Subsample some reads for testing
si -c 16 -t 120
micromamba activate seq_tools

# sample="M1.1-V1"
sample="1226_1"
# assembly_file="mt.assembly.merged.fa.gz"
assembly_file="mgmt.assembly.merged.fa"

sample_dir="/mnt/isilon/projects/ecosystem_biology/MiBiPa/IMP/${sample}/run1"
# sample_dir="/mnt/isilon/projects/ecosystem_biology/user/ohickl/data_backup/must/MUST-totalRNA-plasma/imp3/${sample}"

# test_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/multi_imp_tests/test_imp3_data/must_bloodrna_m1_1_v1_subsample"
test_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/multi_imp_tests/test_imp3_data/mibipa_1226_1_subsample"

# for read_file in ${sample_dir}/Preprocessing/m[g,t].r[1,2].fq; do
#     seqkit sample -p 0.01 -o ${test_dir}/${read_file##*/} ${read_file}
# done

seqkit sample -p 0.1 -o ${test_dir}/${assembly_file} ${sample_dir}/Assembly/${assembly_file}

micromamba activate /mnt/isilon/projects/ecosystem_biology/local_tools/envs/snakemake_env_dev

imp_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"

snakemake --profile ${imp_dir}/cluster_profile_test_aion \
          --snakefile ${imp_dir}/Snakefile  \
          --cluster-cancel scancel  \
          --cluster-status ${imp_dir}/runscripts/status.py  \
          --verbose  \
          --configfile ${imp_dir}/config/config.imp.hybrid_test.yaml -n > ${imp_dir}/imp3_test_0001.log 2>&1
