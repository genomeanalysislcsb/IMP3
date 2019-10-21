#! /bin/bash

#$ -S /bin/bash
#$ -o IMPmgh.stdout
#$ -e IMPmgh.stderr
#$ -l h_rt=24:00:00,h_vmem=6G
#$ -pe smp 3
#$ -binding linear:3
#$ -N IMPmgh.rep


date
module load python/2.7.13-1-ecogen
module load R/3.4.2-1
module load graphviz/2.40.1-1
module load anaconda/3/5.0.1-1
PATH=$PATH:/home/heintzbu/.conda/envs/my_root/bin

CONFIGFILE="report.config.json" /home/heintzbu/.conda/envs/my_root/bin/snakemake -s /data/idiv-biu/heintzbu/TOOLS/IMPscripts/IMP/Snakefile -j 3 report.done
date
