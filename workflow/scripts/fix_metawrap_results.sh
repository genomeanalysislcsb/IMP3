#! /bin/bash -i

# SBATCH --output=/home/$USER/personal/%x-%j.stdout
# SBATCH --error=/home/$USER/personal/%x-%j.stderr
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G

SAMPLE_PATH=$1
cd $SAMPLE_PATH

IMP_PATH=/zfs/omics/projects/metatools/DEVEL/IMP3_versions/220908_IMP3_binning

conda activate ${IMP_PATH}/conda/0adb9bb71a174157a376f82144860244

SOFT=${CONDA_PREFIX}/bin/metawrap-scripts

echo "Re-running CheckM on bins"
mkdir -p Binning/MetaWrap/bins_recheck.tmp

checkm lineage_wf -x fasta Binning/MetaWrap/bins Binning/MetaWrap/bins/checkm -t $SLURM_CPUS_PER_TASK --tmpdir Binning/MetaWrap/bins_recheck.tmp --pplacer_threads 1

if [[ ! -s Binning/MetaWrap/bins/checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
rm -r Binning/MetaWrap/bins_recheck.tmp
        
${SOFT}/summarize_checkm.py Binning/MetaWrap/bins/checkm/storage/bin_stats_ext.tsv  > Binning/MetaWrap/bins.stats
if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi

rm -r Binning/MetaWrap/bins/checkm

conda deactivate

${IMP_PATH}/workflow/scripts/fix_metawrap_results.py Binning/MetaWrap/bins.stats Stats/mg/mg.bins.tsv 

rm Stats/all*

conda activate ${IMP_PATH}/conda/f7b0492ec41a01063bc8974fde389316

Rscript ${IMP_PATH}/workflow/scripts/fix_metawrap_results.R

conda deactivate


