#!/bin/bash

# This script will generate a config file for each sample in the sample list
imp_dir=$1
conf_template=$2
conf_dir=$3
sample_list=$4
sample_list_has_header=${5:-"true"}
imp_run_out_dir=${6:-"./"}

start_line=1

# Check if there's a header to skip
if [ "$sample_list_has_header" == "true" ]; then
  start_line=2
fi

while IFS=$'\t' read -r id mg1 mg2 mg_se mt1 mt2 mt_se; do
  # If the line is just a newline, skip
  if [ -z "$id" ]; then
    continue
  fi

  echo "Printing config for sample ${id}"

  sed_commands="s|__sample_id__|${id}|g"
  sed_commands+=";s|__out_dir__|${imp_run_out_dir}|g"

  # Construct MG input
  mg_input=""
  [ "$mg1" != "NA" ] && mg_input+="$mg1"
  [ "$mg2" != "NA" ] && mg_input+=${mg_input:+\ }"$mg2" # Add space if mg_input is not empty
  [ "$mg_se" != "NA" ] && mg_input+=${mg_input:+\ }"$mg_se" # Add space if mg_input is not empty
  sed_commands+=";s|__input_mg__|${mg_input}|g"

  # Construct MT input
  mt_input=""
  [ "$mt1" != "NA" ] && mt_input+="$mt1"
  [ "$mt2" != "NA" ] && mt_input+=${mt_input:+\ }"$mt2" # Add space if mt_input is not empty
  [ "$mt_se" != "NA" ] && mt_input+=${mt_input:+\ }"$mt_se" # Add space if mt_input is not empty
  sed_commands+=";s|__input_mt__|${mt_input}|g"

  sed -e "${sed_commands}" ${conf_template} > ${conf_dir}/${id}.yaml
  mkdir -p ${imp_run_out_dir}/tmp/${id}

done < <(tail -n +${start_line} ${sample_list})
