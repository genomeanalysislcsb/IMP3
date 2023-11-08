# Create imp config from template
proj_dir="/mnt/lscratch/users/ohickl/MUST-totalRNA-plasma"
imp_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"
conf_template="${imp_dir}/config/config.imp.template.yaml"
conf_dir="${proj_dir}/imp_configs"
sample_list="/mnt/lscratch/users/ohickl/MUST-totalRNA-plasma/samples.tsv"

out_dir="/mnt/lscratch/users/ohickl/MUST-totalRNA-plasma/imp3"

# while IFS=$'\t' read -r id mg1 mg2 mt1 mt2; do
#     -e "s|__input_mt__|${mt1} ${mt2}|g" \
# echo "ID: ${id}; MG_R1: ${mg1}; MG_R2: ${mg2}; MG_SE: ${mg_s}; MT_R1: ${mt1}; MT_R2: ${mt2}; MT_SE: ${mt_s}"
while IFS=$'\t' read -r id mt1 mt2; do
  echo "Printing ${cluster} config for sample ${id}"
  sed -e "s|__sample_id__|${id}|g" \
      -e "s|__out_dir__|${out_dir}|g" \
      -e "s|__input_mg__||g" \
      -e "s|__input_mt__|${mt1} ${mt2}|g" \
      ${conf_template} > ${conf_dir}/${id}.yaml
  mkdir -p ${out_dir}/tmp/${id}
done < <(tail -n +2 ${sample_list})


# bulk imp3
proj_dir="/mnt/lscratch/users/ohickl/MUST-totalRNA-plasma"
imp_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"
conf_template="${imp_dir}/config/config.imp.template.yaml"
conf_dir="${proj_dir}/imp_configs"
conf_array=(${conf_dir}/*.yaml)

micromamba activate /mnt/isilon/projects/ecosystem_biology/local_tools/envs/snakemake_env

 
