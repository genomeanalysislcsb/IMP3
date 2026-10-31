#!/bin/bash -l

#SBATCH -J imp3_star___SAMPLE_____OME__
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c __THREADS__
#SBATCH --clusters __CLUSTER__
#SBATCH --partition __PARTITION__
#SBATCH --time=00-03:00:00
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
mem_pc_byte="__MEMCORE__"


### For manual ###
# sample="1226_1"
# imp_out="/mnt/isilon/projects/ecosystem_biology/MiBiPa/imp_2023"
# out_dir_root="/mnt/lscratch/users/ohickl/mibipa/star/samples"
# threads=128
# ome="mt"
# mem_pc_byte="1750000000"
### For manual ###

out_dir="${out_dir_root}/${sample}_low_thres"

assembly_org="${imp_out}/${sample}/Assembly/mgmt.assembly.merged.fa"
assembly="${out_dir}/mgmt.assembly.merged.fa"

if [[ ! -f ${assembly_org} ]]; then
    assembly_org="${imp_out}/${sample}/Assembly/${ome}.assembly.merged.fa"
    assembly="${out_dir}/${ome}.assembly.merged.fa"
fi

r_1="${imp_out}/${sample}/Preprocessing/${ome}.r1.preprocessed.fq"
r_2="${imp_out}/${sample}/Preprocessing/${ome}.r2.preprocessed.fq"
r_se="${imp_out}/${sample}/Preprocessing/${ome}.se.preprocessed.fq"

if [[ ( ! -f ${r_1} ) || ( ! -f ${r_2} ) ]]; then
    echo "Read files dont exist. Exiting." && exit 1
fi



if [[ ! -d ${out_dir} ]]; then
    mkdir -p ${out_dir}
fi

log_index="${out_dir}/index.log"
log_paired="${out_dir}/${ome}_mapping_paired.log"
log_single="${out_dir}/${ome}_mapping_single.log"


rsync -avPu ${assembly_org} ${out_dir}

paired_star_mem_out="${out_dir}/${ome}.reads.paired."
single_star_mem_out="${out_dir}/${ome}.reads.single."

paired_mapping_sorted="${out_dir}/${ome}.reads.paired.sorted.bam"
single_mapping_sorted="${out_dir}/${ome}.reads.single.sorted.bam"

paired_unmapped_sorted="${out_dir}/${ome}.reads.paired.unmapped.sorted.bam"
single_unmapped_sorted="${out_dir}/${ome}.reads.single.unmapped.sorted.bam"

SAMHEADER="@RG\\tID:${sample}\\tSM:${ome}"
PREFIX=Assembly/{wildcards.type}.reads

micromamba activate seq_tools

# Get star params
r_max_l=$( seqkit stats -T ${r_1} | cut -f 8 | tail -n +2 )
sjdbOverhang=$( echo "${r_max_l}-1" | bc )

assembly_stats=$( seqkit stats -T ${assembly} | tail -n +2 )
a_contigs=$( echo ${assembly_stats} | cut -d ' ' -f 4 )
a_length=$( echo ${assembly_stats} | cut -d ' ' -f 5 )

genomeChrBinNbits=$( python -c "import math; print( round(min(18, math.log2(max(${a_length}/${a_contigs}, ${r_max_l})))) )" )

limitGenomeGenerateRAM=$( echo "${mem_pc_byte}*${threads}" | bc )

# Index
# --sjdbGTFfile /path/to/annotations.gtf \
# --sjdbOverhang ${sjdbOverhang} \
if [[ ( ! -f ${out_dir}/Genome ) || ( ! -f ${out_dir}/genomeParameters.txt ) || ( ! -f ${out_dir}/chrStart.txt ) || ( ! -f ${out_dir}/chrName.txt ) || ( ! -f ${out_dir}/chrNameLength.txt ) || ( ! -f ${out_dir}/chrLength.txt)  ]]; then  # TODO
    STAR --runThreadN ${threads} \
         --runMode genomeGenerate \
         --genomeDir ${out_dir} \
         --limitGenomeGenerateRAM ${limitGenomeGenerateRAM} \
         --genomeFastaFiles ${assembly} > ${log_index} 2>&1 || { echo "Indexing failed. Exiting." && exit 1; }
fi

# Map paired
    #  --outFilterMultimapNmax 100 \ # default 20, no idea how many
    #  --winAnchorMultimapNmax 150 \ # Should to be >= outFilterMultimapNmax 
c

# Get mapped paired
samtools view --threads ${threads} --excl-flags 1540 --bam ${paired_star_mem_out}Aligned.out.sam 2>> ${log_paired} |  \
    samtools sort --threads ${threads} -O bam -o ${paired_mapping_sorted} >> ${log_paired} 2>&1
samtools index ${paired_mapping_sorted} >> ${log_paired} 2>&1

# Get unmapped paired
# samtools view --threads ${threads} --require-flags 4 --excl-flags 1536 --bam ${paired_star_mem_out}Aligned.out.sam 2>> ${log_paired} |  \
#     samtools sort --threads ${threads} -O bam -o ${paired_unmapped_sorted} >> ${log_paired} 2>&1
# samtools index ${paired_unmapped_sorted} >> ${log_paired} 2>&1

# Map single
STAR --runThreadN ${threads} \
     --runMode alignReads \
     --genomeDir ${out_dir} \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 150 \
     --alignIntronMax 1 \
     --outFileNamePrefix ${single_star_mem_out} \
     --outReadsUnmapped Fastx \
     --outFilterScoreMinOverLread 0.33 \
     --outFilterMatchNminOverLread 0.33 \
     --readFilesIn ${r_se} >> ${log_single} 2>&1

# Get mapped single
samtools view --threads ${threads} --excl-flags 1540 --bam ${single_star_mem_out}Aligned.out.sam 2>> ${log_single} |  \
    samtools sort --threads ${threads} -O bam -o ${single_mapping_sorted} >> ${log_single} 2>&1
samtools index ${single_mapping_sorted} >> ${log_single} 2>&1

# Get unmapped single
# samtools view --threads ${threads} --require-flags 4 --excl-flags 1536 --bam ${single_star_mem_out}Aligned.out.sam 2>> ${log_single} |  \
#     samtools sort --threads ${threads} -O bam -o ${single_unmapped_sorted} >> ${log_single} 2>&1
# samtools index ${single_unmapped_sorted} >> ${log_single} 2>&1

rm -rf ${paired_star_mem_out}Aligned.out.sam ${single_star_mem_out}Aligned.out.sam ${out_dir}/Genome ${out_dir}/SA ${out_dir}/SAindex ${assembly}

# samtools merge --threads ${threads} \
#                - \
#                ${paired_mapping_sorted} \
#                ${single_mapping_sorted} | \
#     samtools sort --threads ${threads} -O bam -o ${paired_mapping_sorted/paired/merged}