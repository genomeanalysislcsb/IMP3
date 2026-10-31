
src_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"
imp_out="/mnt/isilon/projects/ecosystem_biology/MiBiPa/imp_2023"

# star_align.sh manual_imp3_bwa.sh
aligner_script="star_align.sh"

if [[ ${aligner_script} == 'star_align.sh' ]]; then
  out_dir_aligner="star"
else
  out_dir_aligner="bwa"
fi

out_dir="/mnt/lscratch/users/ohickl/mibipa/${out_dir_aligner}/samples"

threads=128
# iris bigmem: 26339000000, aion batch: 1750000000
mem_pc_byte="1750000000"
# iris aion
cluster="aion"
# bigmem batch
partition="batch"

declare -a sample_list=()
# echo "input: ${fname} output: ${BASH_REMATCH[1]}";
for fname in ${imp_out}/*_*; do 
    [[ $fname =~ .*/(.*_.*) ]];
    sample_list+=(${BASH_REMATCH[1]})
done

# sample_list=(1226_1)
# sample_list=(1229_14 1233_7 Batch1_3)

# Run failed etc
# declare -a sample_list=()
# # mt.reads.paired.sorted.bam 
# for i in bwa/samples/*; do
#   sample=${i##*/}
#   if [[ ! -d star/samples/${sample} ]]; then
#     sample=${i##*/}
#     sample_list+=(${sample})
#   fi
# done

# Run 99% too short samples (+ 30% for Batch3_1)
# sample_list=(1228_9 1233_12 1228_16)
sample_list=(1228_9)
## 1228_12 needs bigmem
# sample_list=(1228_12)

echo ${#sample_list[@]}
echo ${sample_list[@]}

ome="mt"

for sample in ${sample_list[@]}; do
  echo "Running mapping for sample ${sample} ${ome}."
  sbatch <(sed -e "s|__SAMPLE__|${sample}|g" \
               -e "s|__IMP_OUT__|${imp_out}|g" \
               -e "s|__OUT_DIR__|${out_dir}|g" \
               -e "s|__OME__|${ome}|g" \
               -e "s|__CLUSTER__|${cluster}|g" \
               -e "s|__PARTITION__|${partition}|g" \
               -e "s|__MEMCORE__|${mem_pc_byte}|g" \
               -e "s|__THREADS__|${threads}|g" ${src_dir}/${aligner_script})
  sleep 0.25
done