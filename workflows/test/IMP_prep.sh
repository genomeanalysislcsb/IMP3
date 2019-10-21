#! /bin/bash -l

# -S /bin/bash
# -o IMPmgh.stdout
# -e IMPmgh.stderr
# -l h_rt=24:00:00,h_vmem=6G
# -pe smp 20
# -binding linear:20
# -N IMPmgh.prep

module load lang/Python/3.5.3-intel-2017a
module load lang/R/3.4.0-intel-2017a-X11-20170314-bare

export imp_root=/work/projects/ecosystem_biology/local_tools/IMPPatrick/IMP
export imp_bin=$imp_root/bin
export imp_lib=$imp_root/lib
export DBPATH=$imp_root/databases

export PATH=$PATH:$imp_bin


date
CONFIGFILE="./config.json" $imp_bin/snakemake -s $imp_root/Snakefile -j 8 analysis.done
date


