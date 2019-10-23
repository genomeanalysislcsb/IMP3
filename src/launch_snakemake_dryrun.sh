
CONFIGFILE=$1

conda activate /work/projects/ecosystem_biology/local_tools/minicondaESB
snakemake -j 1 -s /work/projects/ecosystem_biology/local_tools/IMP3/Snakefile --dryrun --configfile $CONFIGFILE

