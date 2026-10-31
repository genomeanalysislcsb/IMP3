
src_dir="/mnt/isilon/projects/ecosystem_biology/local_tools/IMP3_devel"
imp_out="/mnt/isilon/projects/ecosystem_biology/MiBiPa/imp_2023"
out_dir="/mnt/lscratch/users/ohickl/mibipa/featurecounts/samples"

#mg mt#
omes="mg"
align_suffix='bwa'

bam_root="/mnt/lscratch/users/ohickl/mibipa/${align_suffix}/samples"

threads=64

declare -a sample_list=()

# echo "input: ${fname} output: ${BASH_REMATCH[1]}";
for fname in ${imp_out}/*_*; do 
    [[ $fname =~ .*/(.*_.*) ]];
    sample_list+=(${BASH_REMATCH[1]})
done

# sample_list=(1226_1)
# sample_list=(1229_14 1233_7 Batch1_3)

# OOT Batch2_1 1231_16
sample_list=(Batch2_1)


for sample in ${sample_list[@]}; do
  echo "Running featureCounts for sample ${sample}."
  sbatch <(sed -e "s|__SAMPLE__|${sample}|g" \
               -e "s|__IMP_OUT__|${imp_out}|g" \
               -e "s|__OUT_DIR__|${out_dir}|g" \
               -e "s|__OMES__|${omes}|g" \
               -e "s|__BAM_DIR__|${bam_root}/${sample}|g" \
               -e "s|__OUT_SUFFIX__|${align_suffix}|g" \
               -e "s|__THREADS__|${threads}|g" ${src_dir}/manual_imp_feature_counts.sh)
  sleep 0.25
done