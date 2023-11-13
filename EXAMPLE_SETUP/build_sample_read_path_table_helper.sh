#!/bin/bash

################################################################################
##### Function definitions #####################################################
################################################################################
get_read_affix() {
    case $1 in
        ${r1_file_suffix})
            echo "${r1_file_affix}"
            ;;
        ${r2_file_suffix})
            echo "${r2_file_affix}"
            ;;
        ${se_file_suffix})
            echo "${se_file_affix}"
            ;;
        *)
            echo ""
            ;;
    esac
}


# Define output sample list tsv path
out_sample_tsv="/mnt/lscratch/users/ohickl/imp_guide_test/samples.tsv"

# Directory containing read files. Ensure that this points to your directory
# of read files and that files within are structured uniformly.
# In this case they need to be structured as follows:
#   ../samples_dir/sample_name/<ome>.<read_type>
# e.g. ../my_samples/ERR1234567/mg.r1.preprocessed.fq
samples_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/20230330_IMP3/test"

# Define the naming structure for your read files. Modify these variables
# to match your naming convention if different from below.
# The general structure is: <affix><ome><suffix>
# Ensure to include the separator e.g "." in the variable
# You can leave se_file_name empty if you do not have single-end reads
r1_file_suffix=".R1.fq"
r2_file_suffix=".R2.fq"
se_file_suffix=""

# Leave empty if no affix is used
# If there is an affix, ensure to include the separator e.g "." in the variable
r1_file_affix=""
r2_file_affix=""
se_file_affix=""

# Write header to output TSV file
printf "SampleID\tMG_R1\tMG_R2\tMG_SE\tMT_R1\tMT_R2\tMT_SE\n" > ${out_sample_tsv}

# Iterate over each sample directory in the read files directory
for sample_path in ${samples_dir}/*; do
  echo "Processing directory: ${sample_path}"

  # Skip if not a directory
  if [ ! -d "${sample_path}" ]; then
      continue
  fi

  sample_id=${sample_path##*/}
  sample_line="${sample_id}"

  for ome in "mg" "mt"; do
    for read_type in "${r1_file_suffix}" "${r2_file_suffix}" "${se_file_suffix}"; do
      affix=$(get_read_affix ${read_type})
      # Try lowercase ome first
      read_path="${sample_path}/${affix}${ome}${read_type}"
      if [[ ! -f "${read_path}" ]]; then
        # If not found, try uppercase ome
        read_path="${sample_path}/${affix}${ome^^}${read_type}"
      fi

      if [[ -f "${read_path}" ]]; then
        sample_line+=$'\t'"${read_path}"
      else
        sample_line+=$'\tNA'
      fi
    done
  done

  # If sample_line is just sample_id with tabs, then no read files were found and we skip
  if [[ "${sample_line}" == "${sample_id}"$'\tNA\tNA\tNA\tNA\tNA\tNA' ]]; then
      echo "No read files found for sample: ${sample_id}"
      continue
  fi

  # Append the constructed line to the output TSV file
  printf "%s\n" "${sample_line}" >> ${out_sample_tsv}
done
