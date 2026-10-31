proj_dir="/mnt/lscratch/users/ohickl/MUST-totalRNA-plasma"
# KRAKEN_DB="/mnt/isilon/projects/ecosystem_biology/user/ohickl/imp3_test_db/kraken_struo2_gtdb_207"
KRAKEN_DB="/mnt/isilon/projects/ecosystem_biology/databases/imp3/20230502_kraken_nt"
READ_LEN=50
# CLASSIFICATION_LEVEL="S" # D,P,C,O,F,G,S,S1,etc
# THRESHOLD=3
THREADS=128
KMER_LEN=35
#( mg mt )
omes=( mt )
lvls2run=( D P C O F G S )

out_dir="${proj_dir}/Bracken_reads_rlen${READ_LEN}"

micromamba activate seq_tools

####
# check db
kraken2-inspect --db ${KRAKEN_DB} \
                --threads ${THREADS} \
                --skip-counts

#build db for readlen x
bracken-build -d ${KRAKEN_DB} \
              -t ${THREADS} \
              -k ${KMER_LEN} \
              -l ${READ_LEN}
####

# for SAMPLE in ${proj_dir}/imp_2023/*_*; do
#   echo $sample
# done | less

# SAMPLE="${proj_dir}/imp_2023/Man14122018_2"
# SAMPLE_PATH="${proj_dir}/imp_2023/Man14122018_2"


# Kraken2
conf=0.51
kraken2 --db ${KRAKEN_DB} \
        --confidence ${conf} \
        --threads ${THREADS} \
        --paired \
        --use-names \
        --output ${out_dir}/${SAMPLE}/mt.kraken.output \
        --report {output[1]} \
        {input} > {log} 2>&1

# Contigs makes no sense
# for SAMPLE_PATH in ${proj_dir}/imp3/M*; do
# #   SAMPLE=${SAMPLE_PATH%/}
#   SAMPLE=${SAMPLE_PATH##*/}
#   mkdir -p ${proj_dir}/Bracken_contigs/${SAMPLE}
#   for CLASSIFICATION_LEVEL in ${lvls2run[@]}; do
#     for THRESHOLD in 1 5; do
#       bracken -d ${KRAKEN_DB} \
#               -i ${SAMPLE_PATH}/Analysis/taxonomy/kraken/mgmt.kraken.contig.report \
#               -o ${proj_dir}/Bracken_contigs/${SAMPLE}/${SAMPLE}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.bracken \
#               -r ${READ_LEN} \
#               -l ${CLASSIFICATION_LEVEL} \
#               -t ${THRESHOLD} > ${proj_dir}/Bracken_contigs/logs/${SAMPLE}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.log 2>&1
#     done
#   done
# done | less -RS

# Reads from imp out
mkdir -p ${out_dir}/logs
for SAMPLE_PATH in ${proj_dir}/imp3/M*; do
#   SAMPLE=${SAMPLE_PATH%/}
  SAMPLE=${SAMPLE_PATH##*/}
  mkdir -p ${out_dir}/${SAMPLE}
  for ome in ${omes[@]}; do
    for CLASSIFICATION_LEVEL in ${lvls2run[@]}; do
      for THRESHOLD in 1 5 10; do
        bracken -d ${KRAKEN_DB} \
                -i ${SAMPLE_PATH}/Analysis/taxonomy/kraken/${ome}.kraken.report \
                -o ${out_dir}/${SAMPLE}/${SAMPLE}_${ome}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.bracken \
                -r ${READ_LEN} \
                -l ${CLASSIFICATION_LEVEL} \
                -t ${THRESHOLD} > ${out_dir}/logs/${SAMPLE}_${ome}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.log 2>&1
      done
    done
  done
done

# Reads from manual out
db="struo2_gtdb_207"  # 20230502_kraken_nt struo2_gtdb_207

if [[ ${db} == '20230502_kraken_nt' ]]; then
  KRAKEN_DB="/mnt/isilon/projects/ecosystem_biology/databases/imp3/20230502_kraken_nt"
  READ_LEN=50
else
  KRAKEN_DB="/mnt/isilon/projects/ecosystem_biology/user/ohickl/imp3_test_db/kraken_struo2_gtdb_207"
  READ_LEN=100
fi

echo ${KRAKEN_DB}

conf="0.00"
kraken_outs=(/mnt/isilon/projects/ecosystem_biology/MUST/MUST-plasma/Kraken2/Kraken_0250_conf/M*_mt_reads_${db}_conf_${conf}.kraken.report)

echo ${#kraken_outs[@]}
echo ${kraken_outs[@]}

out_dir="${proj_dir}/Bracken_reads_rlen${READ_LEN}_db_${db}_conf_${conf}"
mkdir -p ${out_dir}/logs

run_bracken () {
  local SAMPLE_PATH=$1
  local KRAKEN_DB=$2
  local READ_LEN=$3
  local omes=( $4 )
  local out_dir=$5
  local SAMPLE=${SAMPLE_PATH%_mt_reads_*}
  local SAMPLE=${SAMPLE##*/}
  mkdir -p ${out_dir}/${SAMPLE}
  for ome in ${omes[@]}; do
    for CLASSIFICATION_LEVEL in D P C O F G S; do
      for THRESHOLD in 1 5 10; do
        bracken -d ${KRAKEN_DB} \
                -i ${SAMPLE_PATH} \
                -o ${out_dir}/${SAMPLE}/${SAMPLE}_${ome}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.bracken \
                -r ${READ_LEN} \
                -l ${CLASSIFICATION_LEVEL} \
                -t ${THRESHOLD} > ${out_dir}/logs/${SAMPLE}_${ome}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.log 2>&1
      done
    done
  done
}

export -f run_bracken


parallel -j 62 \
         "run_bracken {} ${KRAKEN_DB} ${READ_LEN} ${omes[@]} ${out_dir}" ::: ${kraken_outs[@]}

# for SAMPLE_PATH in ${kraken_outs[@]}; do
#   SAMPLE=${SAMPLE_PATH%_mt_reads_*}
#   SAMPLE=${SAMPLE##*/}
#   mkdir -p ${out_dir}/${SAMPLE}
#   for ome in ${omes[@]}; do
#     for CLASSIFICATION_LEVEL in ${lvls2run[@]}; do
#       for THRESHOLD in 1 5 10; do
#         bracken -d ${KRAKEN_DB} \
#                 -i ${SAMPLE_PATH} \
#                 -o ${out_dir}/${SAMPLE}/${SAMPLE}_${ome}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.bracken \
#                 -r ${READ_LEN} \
#                 -l ${CLASSIFICATION_LEVEL} \
#                 -t ${THRESHOLD} > ${out_dir}/logs/${SAMPLE}_${ome}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.log 2>&1
#       done
#     done
#   done
# done

# combine reports
db="20230502_kraken_nt"  # 20230502_kraken_nt struo2_gtdb_207

if [[ ${db} == '20230502_kraken_nt' ]]; then
  READ_LEN=50
else
  READ_LEN=100
fi

conf="0.00"

CLASSIFICATION_LEVEL="S" # D,P,C,O,F,G,S,S1,etc
lvls2run=( D P C O F G S )

# out_dir="${proj_dir}/Bracken_reads_rlen${READ_LEN}"
out_dir="${proj_dir}/Bracken_reads_rlen${READ_LEN}_db_${db}_conf_${conf}"


CLASSIFICATION_LEVEL=G
THRESHOLD=10
ome=mt

for CLASSIFICATION_LEVEL in ${lvls2run[@]}; do
  for THRESHOLD in 1 5 10; do
    reports=(${out_dir}/*/*_${ome}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.bracken)
    names=(${reports[@]/${out_dir}\//})
    names=(${names[@]/\/*_${ome}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.bracken/})
    printf -v names_csv '%s,' "${names[@]}"
    names_csv="${names_csv%,}"

    combine_bracken_outputs.py --files ${reports[@]} \
                              --names ${names_csv[@]} \
                              --output ${out_dir}/all_samples_${ome}_thresh_${THRESHOLD}_level_${CLASSIFICATION_LEVEL}.bracken
  done
done

