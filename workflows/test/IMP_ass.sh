#! /bin/bash

#$ -S /bin/bash
#$ -o IMPmgh.stdout
#$ -e IMPmgh.stderr
#$ -l h_rt=48:00:00,h_vmem=18G,highmem
#$ -pe smp 2
#$ -binding linear:2
#$ -N IMPmgh.ass


date
module load bedtools/2.20.1-1
module load samtools/1.1-2
module load bwa/0.7.12-1
module load anaconda/3/5.0.1-1

CONFIGFILE="assembly.config.json" /home/heintzbu/.conda/envs/my_root/bin/snakemake -s /data/idiv-biu/heintzbu/TOOLS/IMPscripts/IMP/Snakefile -j 2 assembly.done
date
