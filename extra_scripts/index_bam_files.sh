
micromamba activate seq_tools

parallel -j 128 --dry-run \
         'samtools index {}' ::: */samples/*/*.bam