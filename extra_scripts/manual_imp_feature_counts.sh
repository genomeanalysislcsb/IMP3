#!/bin/bash -l

#SBATCH -J imp3_featurecounts___SAMPLE__
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c __THREADS__
#SBATCH --time=00-03:00:00
#SBATCH -o "__OUT_DIR__/../slurm_logs/%x_%j.out"

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"


sample="__SAMPLE__"
imp_out="__IMP_OUT__"
out_dir="__OUT_DIR__"
threads=__THREADS__
omes="__OMES__"
bam_dir="__BAM_DIR__"
out_suffix="__OUT_SUFFIX__"

### For manual ###
# sample="1226_1"
# imp_out="/mnt/isilon/projects/ecosystem_biology/MiBiPa/imp_2023"
# out_dir="/mnt/lscratch/users/ohickl/mibipa/imp_fc_fraction/samples"
# threads=64
# omes="mt"
# bam_dir="/mnt/lscratch/users/ohickl/mibipa/star/samples/1226_1"
# out_suffix="test"
### For manual ###

omes_arr=( ${omes} )

imp_sample_out="${imp_out}/${sample}"
out_sample_dir="${out_dir}/${sample}_${out_suffix}"


mkdir -p ${out_sample_dir}/logs ${out_sample_dir}/gffs

fc_env="/mnt/isilon/projects/ecosystem_biology/local_tools/envs/imp3/4ded0c70e43c7c49116deea8197728af_"

micromamba activate ${fc_env}

dbs=(arcog cazy cog eggnog enzyme_ec go is_essential_gene kegg_brite kegg_ko kegg_module kegg_pathway kegg_rclass kegg_reaction pfam tcdb tigrfam)
features=(CDS rRNA)

gff="${imp_sample_out}/Analysis/annotation/annotation_CDS_RNA_hmms.gff"

int_dir="${imp_sample_out}/Analysis/annotation/intermediary.tar.gz"

rsync -avPu ${int_dir} ${out_sample_dir} \
  || rsync -avPu --include="consensus_annotation.mantis.*.gff" ${int_dir%%.tar.gz}/mantis/* ${out_sample_dir}/gffs \
  || exit "Error: Could not get imp3 annotation output."

if [[ -f ${out_sample_dir}/intermediary.tar.gz ]]; then
    tar -xvf ${out_sample_dir}/intermediary.tar.gz \
        --directory ${out_sample_dir}/gffs \
        --wildcards "Analysis/annotation/intermediary/mantis/consensus_annotation.mantis.*.gff" \
        --strip-components 4 \
      && rm -rf ${out_sample_dir}/intermediary.tar.gz
fi

if ! test -n "$(find ${out_sample_dir}/gffs -maxdepth 1 -name 'consensus_annotation.mantis.*.gff' -print -quit)"; then
  exit "Error: Could not find copied imp3 mantis annotation output." 
fi

for ome in ${omes[@]}; do
    if [[ ${ome} == 'mg' ]]; then
        stranding=0,0
    else
        stranding="2,0"
    fi

    # mapping="${imp_sample_out}/Assembly/${ome}.reads.sorted.bam"

    mapping_pair="${bam_dir}/${ome}.reads.paired.sorted.bam"
    mapping_single="${bam_dir}/${ome}.reads.single.sorted.bam"
    mapping=${mapping_pair}
    mappings="${mapping_pair} ${mapping_single}"

    # mapping="${test_bam_dir}/${ome}.reads.merged.sorted.bam"
    # mappings=${mapping}

    if [[ -f ${mapping} ]]; then
        for flag_setup in '-p -f -O --fraction -M --fraction' '-p -O --fraction -M --fraction'; do
            if [[ ${flag_setup} == '-p -O --fraction -M --fraction' ]]; then
                suffix="feature_count"
                fc_out_dir="${out_sample_dir}"
            else
                suffix="feature_count_per_contig"
                fc_out_dir="${out_sample_dir}/${ome}"
            fi

            if [[ ! -d ${fc_out_dir} ]]; then
                mkdir -p ${fc_out_dir}
            fi

            for db in ${dbs[@]}; do
                if [[ ! -f ${fc_out_dir}/${ome}.mantis.${db}.counts.tsv ]]; then
                    echo "Sample ${sample}: Running featureCounts with flags '${flag_setup}' on ${ome} with db ${db}."
                    featureCounts ${flag_setup} \
                                -t CDS \
                                -g ${db} \
                                -o ${fc_out_dir}/${ome}.mantis.${db}.counts.tsv \
                                -s ${stranding} \
                                -a ${out_sample_dir}/gffs/consensus_annotation.mantis.${db}.gff \
                                -T ${threads} \
                                ${mappings} > ${out_sample_dir}/logs/stats_${suffix}.${ome}.mantis.${db}.log 2>&1 \
                    || touch ${fc_out_dir}/${ome}.mantis.${db}.counts.tsv ${fc_out_dir}/analysis_${suffix}.${ome}.${db}.impossible
                else
                    echo "${fc_out_dir}/${ome}.mantis.${db}.counts.tsv already exists. Skipping."
                fi
            done
        
            for feature in ${features[@]}; do
                if [[ ! -f ${fc_out_dir}/${ome}.${feature}_counts.tsv ]]; then
                    echo "Sample ${sample}: Running featureCounts with flag '${flag_setup}' on ${ome} with feature ${feature}."
                    featureCounts ${flag_setup} \
                                -t ${feature} \
                                -g ID \
                                -o ${fc_out_dir}/${ome}.${feature}_counts.tsv \
                                -s ${stranding} \
                                -a ${gff} \
                                -T ${threads} \
                                ${mappings} > ${out_sample_dir}/logs/analysis_${suffix}_gene.${feature}.${ome}.genes.log 2>&1 \
                    || touch ${fc_out_dir}/${ome}.${feature}_counts.tsv ${fc_out_dir}/analysis_${suffix}_gene.${ome}.${feature}.impossible
                else
                    echo "${fc_out_dir}/${ome}.${feature}_counts.tsv already exists. Skipping."
                fi
            done
        done
    fi
done
