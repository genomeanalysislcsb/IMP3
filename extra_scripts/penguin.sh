# Activate conda env
micromamba activate seq_tools

###### VARS ##############################################################################
threads=128

min_id=0.99  
min_cov=0.90
cluster_mode=2
coverage_mode=1

scratch_dir="/mnt/lscratch/users/ohickl"

read_dir="/mnt/lscratch/users/ohickl/tests/imp_test/Preprocessing"
r1="${read_dir}/mt.r1.trimmed.rrna.fq.gz"
r2="${read_dir}/mt.r2.trimmed.rrna.fq.gz"
se="${read_dir}/mt.se.trimmed.rrna.fq.gz"
# r1="${read_dir}/mt.r1.trimmed.rrna.fq.gz"
# r2="${read_dir}/mt.r2.trimmed.rrna.fq.gz"


proj_dir="${scratch_dir}/assembly_tests/penguin"
out_dir="${proj_dir}/penguin_rna_001"


out_prefix="${out_dir}/clustered"
assembly_fasta="${out_dir}/assembly.fasta"
tmp_dir="/dev/shm"

assembly_fasta_tmp="${tmp_dir}/assembly.fasta"

log="${out_dir}/assembly.log"

# Ensure output directory exists
mkdir -p ${out_dir}

# # To use se reads we first need to merge paired reads since penguin does not support mixing paired and single end reads
# # Merge paired reads
# flash2 ${r1} \
#        ${r2} \
#        --interleaved-output \
#        --output-prefix merged \
#        --output-directory ${tmp_dir} \
#        --threads ${threads} > ${log} 2>&1

###### Assembly ##########################################################################
echo "###### Running penguin on paired reads ######" > ${log}
micromamba run -n plass \
    penguin guided_nuclassemble \
        ${r1} \
        ${r2} \
        ${assembly_fasta_tmp/.fasta/}_paired.fasta \
        ${tmp_dir} \
        --min-length 21 \
        --use-all-table-starts 1 \
        --min-contig-len 200 \
        --contig-output-mode 0 \
        --num-iterations aa:5,nucl:15 \
        --threads 10 \
        --max-seq-len 9999999 >> ${log} 2>&1

echo "###### Running penguin on single end reads with paired contigs ######" >> ${log}
# First we make a fasta out of the se fq file
seqkit fq2fa ${se} > ${assembly_fasta_tmp/.fasta/}_se_reads.fasta 2>> ${log}
# Then we concat with the paired assembly
cat ${assembly_fasta_tmp/.fasta/}_paired.fasta \
    ${assembly_fasta_tmp/.fasta/}_se_reads.fasta > ${assembly_fasta_tmp/.fasta/}_se.fasta 2>> ${log} \
    && rm ${assembly_fasta_tmp/.fasta/}_se_reads.fasta >> ${log} 2>&1

micromamba run -n plass \
    penguin guided_nuclassemble \
        ${assembly_fasta_tmp/.fasta/}_se.fasta \
        ${assembly_fasta_tmp} \
        ${tmp_dir} \
        --min-length 21 \
        --use-all-table-starts 1 \
        --min-contig-len 200 \
        --contig-output-mode 0 \
        --num-iterations aa:5,nucl:15 \
        --threads 10 \
        --max-seq-len 9999999 >> ${log} 2>&1

# test iteration, assembly fasta 10 times, each time taking the output and using it as input
# for i in {1..10}; do
#     next_iteration=$((${i}+1))
#     micromamba run -n plass \
#         penguin guided_nuclassemble \
#             ${assembly_fasta_tmp}${i} \
#             --min-length 21 \
#             --use-all-table-starts 1 \
#             --min-contig-len 200 \
#             --contig-output-mode 0 \
#             --num-iterations aa:5,nucl:15 \
#             --threads 10 \
#             --max-seq-len 9999999 \
#             ${assembly_fasta_tmp}${next_iteration} \
#             ${tmp_dir}
# done


# # Combine paired and single end assemblies
# cat ${assembly_fasta_tmp/.fasta/}_paired.fasta \
#     ${assembly_fasta_tmp/.fasta/}_se.fasta > ${assembly_fasta_tmp} >2 ${log}

###### Clustering ########################################################################
# Cluster at 99% identity as in paper
echo "###### Clustering at ${min_id} identity and ${min_cov} coverage ######" >> ${log}
micromamba run -n mmseqs \
    mmseqs easy-cluster \
        ${assembly_fasta_tmp} \
        ${out_prefix} \
        ${tmp_dir} \
        --min-seq-id ${min_id} \
        --cov-mode ${coverage_mode} \
        -c ${min_cov} \
        --cluster-mode ${cluster_mode} \
        --threads ${threads} >> ${log} 2>&1

# Compress large fasta
pigz -p ${threads} ${out_prefix}_all_seqs.fasta >> ${log} 2>&1

# Move everything to final location
rsync -av ${tmp_dir}/ ${out_dir}/ >> ${log} 2>&1

# Cleanup
rm -rf ${tmp_dir}/* >> ${log} 2>&1

