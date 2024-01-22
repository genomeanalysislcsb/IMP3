#!/bin/bash

samples_path="/work/projects/ecosystem_biology/MUST/MUST-totalRNA-plasma/Raw/M*1.fq"

sample_tsv="/mnt/lscratch/users/ohickl/MUST-totalRNA-plasma/samples.tsv"
# printf "SampleID\tMG_R1\tMG_R2\tMG_SE\tMT_R1\tMT_R2\tMT_SE\n" > ${sample_tsv}
printf "SampleID\tMT_R1\tMT_R2\n" > ${sample_tsv}

#   se=$(python -c "print(\"${r1}\".replace('.1.fq', '.se.fq'))")
for r1 in ${samples_path}; do
  sample_id=$(python -c "print(\"${r1}\".split('/')[-1].replace('.1.fq', ''))")
  r2=$(python -c "print(\"${r1}\".replace('.1.fq', '.2.fq'))")
  sample_line="${sample_id}\t${r1}\t${r2}\n"
  printf ${sample_line} >> ${sample_tsv}
done