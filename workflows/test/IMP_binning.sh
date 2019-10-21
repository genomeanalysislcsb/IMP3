#! /bin/bash

#$ -S /bin/bash
#$ -o IMPmgh.stdout
#$ -e IMPmgh.stderr
#$ -l h_rt=24:00:00,h_vmem=6G
#$ -pe smp 6
#$ -binding linear:6
#$ -N IMPmgh.bin


date
module load python/2.7.13-1-ecogen
module load R/3.4.2-1
module load hmmer/3.1b1-1
module load pcre/8.38-1
module load perl/5.22/0-1
module load anaconda/3/5.0.1-1

CONFIGFILE="binning.config.json" /home/heintzbu/.conda/envs/my_root/bin/snakemake -s /data/idiv-biu/heintzbu/TOOLS/IMPscripts/IMP/Snakefile -j 6 binning.done
date
