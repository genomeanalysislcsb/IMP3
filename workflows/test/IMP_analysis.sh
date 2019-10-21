#! /bin/bash

#$ -S /bin/bash
#$ -o IMPmgh.stdout
#$ -e IMPmgh.stderr
#$ -l h_rt=48:00:00,h_vmem=6G
#$ -pe smp 20
#$ -binding linear:20
#$ -N IMPmgh.ana


date
module load bedtools/2.20.1-1
module load vcftools/0.1.12-2
module load bcftools/1.1-1
module load samtools/0.1.19-1
module load htslib/1.3-1
module load python/2.7.13-1-ecogen
module load R/3.4.2-1
module load fastqc/0.11.2-1
module load perl/5/22/0-1
module load hmmer/3.1b1-1
module load anaconda/3/5.0.1-1

CONFIGFILE="analysis.config.json" /home/heintzbu/.conda/envs/my_root/bin/snakemake -s /data/idiv-biu/heintzbu/TOOLS/IMPscripts/IMP/Snakefile -j 20 analysis.done
date
