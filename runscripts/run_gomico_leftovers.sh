gomico_out_path="/mnt/isilon/projects/ecosystem_biology/GoMiCo"
imp_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"
conf_dir="${gomico_out_path}/sample_imp_configs_v2"
conf_template="${imp_dir}/config/config.imp.template.gomico.yaml"


# Get sample array
declare -a sample_list=()
for sample_path in ${gomico_out_path}/G*; do 
    [[ ${sample_path} =~ .*/(G.*) ]];
    sample_list+=(${BASH_REMATCH[1]})
done

# Write sample imp configs
for sample in ${sample_list[@]}; do
  echo "Printing imp config for sample ${sample}"
  sed -e "s|__sample_id__|${sample}|g" \
      ${conf_template} > ${conf_dir}/${sample}.yaml
done


# Get config array
declare -a conf_array=()
for sample_path in ${gomico_out_path}/G*; do 
    if [[ ! -f ${sample_path}/Stats/all_bin_stats.tsv ]]; then
      [[ ${sample_path} =~ .*/(G.*) ]];
      conf_array+=(${conf_dir}/${BASH_REMATCH[1]}.yaml)
    fi
done

echo ${conf_array[@]}
echo ${#conf_array[@]}

# Run imp
parallel -j 10 \
         "snakemake --profile ${imp_dir}/cluster_profile_gomico " \
                   "--snakefile ${imp_dir}/Snakefile " \
                   "--cluster-status ${imp_dir}/runscripts/status.py " \
                   "--configfile {} --unlock; " \
         "snakemake --profile ${imp_dir}/cluster_profile_gomico " \
                   "--snakefile ${imp_dir}/Snakefile " \
                   "--cluster-status ${imp_dir}/runscripts/status.py " \
                   "--verbose " \
                   "--configfile {} > ${gomico_out_path}/logs/{/.}_hybrid_01.log 2>&1" ::: "${conf_array[@]}" > ${gomico_out_path}/logs/run_gomico_leftovers_01.log 2>&1