
CONFIGFILE=$1

conda activate /work/projects/ecosystem_biology/local_tools/minicondaESB

snakemake -j 20 -s /work/projects/ecosystem_biology/local_tools/IMP3/Snakefile --cluster "sbatch --time={params.runtime} {params.mem} {params.partition} -n {threads}" --configfile $CONFIGFILE --use-conda --conda-prefix /work/projects/ecosystem_biology/local_tools/IMP3/conda

