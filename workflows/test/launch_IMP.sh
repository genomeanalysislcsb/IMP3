#! /bin/bash

oarsub -lcore=8,walltime=96 -t bigmem --project project_biocore -p "network_address='gaia-74'" -n testIMP  ./IMP_prep.sh
#qsub -cwd -n testIMP -hold_jid IMPmgh.prep ./IMP_ass.sh
#qsub -cwd -n testIMP -hold_jjd IMPmgh.ass ./IMP_analysis.sh
#qsub -cwd -n testIMP -hold_jid IMPmgh.ana ./IMP_binning.sh
#qsub -cwd -n testIMP -hold_jid IMPmgh.bin ./IMP_report.sh

