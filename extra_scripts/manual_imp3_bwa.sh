#!/bin/bash -l

#SBATCH -J imp3_bwa-mem2___SAMPLE_____OME__
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c __THREADS__
#SBATCH --clusters __CLUSTER__
#SBATCH --partition __PARTITION__
#SBATCH --time=00-04:00:00
#SBATCH -o "__OUT_DIR__/../slurm_logs/%x_%j.out"

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"


sample="__SAMPLE__"
imp_out="__IMP_OUT__"
out_dir_root="__OUT_DIR__"
threads=__THREADS__
ome="__OME__"

### For manual ###
# sample="1226_1"
# imp_out="/mnt/isilon/projects/ecosystem_biology/MiBiPa/imp_2023"
# out_dir_root="/mnt/lscratch/users/ohickl/mibipa/bwa/samples"
# threads=64
# ome="mt"
### For manual ###

assembly_org="${imp_out}/${sample}/Assembly/mgmt.assembly.merged.fa"

if [[ ! -f ${assembly_org} ]]; then
    assembly_org="${imp_out}/${sample}/Assembly/${ome}.assembly.merged.fa"
fi

r_1="${imp_out}/${sample}/Preprocessing/${ome}.r1.preprocessed.fq"
r_2="${imp_out}/${sample}/Preprocessing/${ome}.r2.preprocessed.fq"
r_se="${imp_out}/${sample}/Preprocessing/${ome}.se.preprocessed.fq"

if [[ ( ! -f ${r_1} ) || ( ! -f ${r_2} ) ]]; then
    echo "Read files dont exist. Exiting." && exit 1
fi

out_dir="${out_dir_root}/${sample}"

if [[ ! -d ${out_dir} ]]; then
    mkdir -p ${out_dir}
fi

log_index="${out_dir}/index.log"
log_paired="${out_dir}/${ome}_mapping_paired.log"
log_single="${out_dir}/${ome}_mapping_single.log"


rsync -avPu ${assembly_org} ${out_dir}

assembly="${out_dir}/mgmt.assembly.merged.fa"

paired_bwa_mem_out="${out_dir}/${ome}.reads.paired.sam"
single_bwa_mem_out="${out_dir}/${ome}.reads.single.sam"

paired_mapping_sorted="${out_dir}/${ome}.reads.paired.sorted.bam"
single_mapping_sorted="${out_dir}/${ome}.reads.single.sorted.bam"

paired_unmapped_sorted="${out_dir}/${ome}.reads.paired.unmapped.sorted.bam"
single_unmapped_sorted="${out_dir}/${ome}.reads.single.unmapped.sorted.bam"

SAMHEADER="@RG\\tID:${sample}\\tSM:${ome}"
PREFIX=Assembly/{wildcards.type}.reads

micromamba activate seq_tools

# Index
if [[ ( ! -f ${assembly}.0123 ) || ( ! -f ${assembly}.amb ) || ( ! -f ${assembly}.ann ) || ( ! -f ${assembly}.bwt.2bit.64 ) || ( ! -f ${assembly}.pac ) ]]; then
    bwa-mem2 index ${assembly} > ${log_index} 2>&1 
fi

# Map paired
bwa-mem2 mem -v 1 -t ${threads} -M -R ${SAMHEADER} -a ${assembly} ${r_1} ${r_2} 2>> ${log_paired} 1>> ${paired_bwa_mem_out}

# Get mapped paired
samtools view --threads ${threads} --excl-flags 1540 --bam ${paired_bwa_mem_out} 2>> ${log_paired} |  \
    samtools sort --threads ${threads} -O bam -o ${paired_mapping_sorted} >> ${log_paired} 2>&1 
samtools index ${paired_mapping_sorted} >> ${log_paired} 2>&1

# Get unmapped paired
samtools view --threads ${threads} --unmap --require-flags 4 --excl-flags 1536 --bam ${paired_bwa_mem_out} 2>> ${log_paired} |  \
    samtools sort --threads ${threads} -O bam -o ${paired_unmapped_sorted} >> ${log_paired} 2>&1 
samtools index ${paired_unmapped_sorted} >> ${log_paired} 2>&1

# Map single
bwa-mem2 mem -v 1 -t ${threads} -M -R ${SAMHEADER} -a ${assembly} ${r_se} 2>> ${log_single} 1>> ${single_bwa_mem_out}

# Get mapped single
samtools view --threads ${threads} --excl-flags 1540 --bam ${single_bwa_mem_out} 2>> ${log_single} |  \
    samtools sort --threads ${threads} -O bam -o ${single_mapping_sorted} >> ${log_single} 2>&1 
samtools index ${single_mapping_sorted} >> ${log_single} 2>&1

# Get unmapped single
samtools view --threads ${threads} --unmap --require-flags 4 --excl-flags 1536 --bam ${single_bwa_mem_out} 2>> ${log_single} |  \
    samtools sort --threads ${threads} -O bam -o ${single_unmapped_sorted} >> ${log_single} 2>&1 
samtools index ${single_unmapped_sorted} >> ${log_single} 2>&1

rm -rf ${paired_bwa_mem_out} ${single_bwa_mem_out} ${assembly}

# samtools merge --threads ${threads} \
#                - \
#                ${paired_mapping_sorted} \
#                ${single_mapping_sorted} | \
#     samtools sort --threads ${threads} -O bam -o ${paired_mapping_sorted/paired/merged}