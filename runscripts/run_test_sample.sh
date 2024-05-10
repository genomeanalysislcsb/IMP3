#!/bin/bash -l

# Run a test sample with IMP
env_manager="micromamba" # conda mamba micromamba

SM_ENV_DIR="/mnt/isilon/projects/ecosystem_biology/local_tools/envs/snakemake_env_dev"
${env_manager} activate ${SM_ENV_DIR} || { echo "Error: Failed to activate env ${SM_ENV_DIR}"; exit 1; }

imp_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"
snakefile="${imp_dir}/Snakefile"

proj_path="/mnt/lscratch/users/ohickl/prospectomics/tests"
sample="imp_test"

s_dir="${proj_path}/${sample}"
prep_path="${s_dir}/Preprocessing"

# Get current date and time
datetime=$(date '+%Y-%m-%d_%H-%M-%S')

log_path="${proj_path}/imp3_test_${datetime}.log"

dryrun=""


# Unlock to ensure no issues with snakemake
snakemake --snakefile ${snakefile} \
          --unlock \
          --quiet all

# Ensure preprocessing is finished
omes=("mg" "mt")
reads=("r1" "r2" "se")
files_to_make=("fq" "trimmed.fq" "trimmed.rrna.fq" "trimmed__rna_filtered.fq" "trimmed__rrna_proc.fq" "trimmed__rrna_proc__PhiX_filtered.fq" "trimmed__rrna_proc__PhiX_filtered__T2TCHM13v2genomic_filtered.fq" "trimmed__rrna_proc__PhiX_filtered__T2TCHM13v2genomic_filtered__T2TCHM13v2transcripts_filtered.fq")

if [[ -d "${prep_path}" ]]; then
  # Check for non-empty '*preprocessed*' files. If none are found, delete empty files.
  if ! find "${prep_path}" -type f -name '*preprocessed*' -size +0c | read; then
      echo "No non-empty '*preprocessed*' files found. Deleting empty files..." > ${log_path}
      # Find and delete empty files in prep_path
      find "${prep_path}" -type f -empty -delete >> ${log_file} 2>&1
  else
      echo "Non-empty '*preprocessed*' files found. Touching temp files..." >> ${log_path}
    for ome in "${omes[@]}"; do
      for read in "${reads[@]}"; do
        for file in "${files_to_make[@]}"; do
          echo "Creating ${prep_path}/${ome}.${read}.${file} ..." >> ${log_path}
          touch --date="2024-03-15 09:30:00" "${prep_path}/${ome}.${read}.${file}" >> ${log_path} 2>&1
        done
      
        if [ -f "${prep_path}/${ome}.${read}.preprocessed.fq" ]; then
          echo "Updating date for ${prep_path}/${ome}.${read}.preprocessed.fq ..." >> ${log_path}
          touch -m --date="2024-03-15 09:31:00" "${prep_path}/${ome}.${read}.preprocessed.fq" >> ${log_path} 2>&1
        fi
      done
    done
  fi
fi

snakemake --profile ${imp_dir}/cluster_profile_test \
          --snakefile ${snakefile} \
          --verbose  \
          --until comebin_bin \
          ${dryrun} >> ${log_path} 2>&1
