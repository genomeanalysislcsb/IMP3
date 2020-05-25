# General

This branch has been derived from the branch `IMPiris`.

- [documentation](https://imp3.readthedocs.io/en/impiris/) for the branch `IMPiris`

# Set-up

**Notes**
- `conda` environment: copy of `submodules/PathoFact/envs/PathoFact.yaml`
    - file: `config/conda.imp.yaml`
- database files: mostly copied from `/work/projects/ecosystem_biology/local_tools/IMP3/databases/`


```bash
# clone project
git clone -b IMPiris_VG ssh://git@git-r3lab-server.uni.lu:8022/IMP/imp3.git IMP3_VG
# update submodules
git submodule update --init --recursive

# create conda env
conda create env -f config/conda.imp.yaml

# databases
DB_PATH="/work/projects/ecosystem_biology/local_tools/IMP3/databases"
mkdir -p databases && cd databases
# links
find ${DB_PATH} -maxdepth 1 -mindepth 1 -type l | while read f; do ln -s $(realpath ${f}) $(basename ${f}); done
# files
find ${DB_PATH} -maxdepth 1 -mindepth 1 -type f -exec rsync {} . \;
# dirs
find ${DB_PATH} -maxdepth 1 -mindepth 1 -type d -exec rsync -avP {} . \;
# TIGRPFAM
cd hmm && ln -s /work/projects/ecosystem_biology/data/hmm/TIGRPFAM TIGRPFAM
```

Follow installation instructions for `PathoFact`: see the `README.md` in `submodules/PathoFact/`.

# Testing the pipeline

Using dada from `test/`

```bash
# dry-run (-n)
conda activate IMP_PF
snakemake -s Snakefile -rp --cores 1 --configfile test/config.imp.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/IMP -n
```

Example for `test/config.imp.yaml`:
```yaml
steps: "preprocessing assembly analysis summary"
preprocessing_filtering: true
tmp_dir: tmp
raws: 
  Metagenomics: "PathToIMP3/test/Reads/test.mg.1.fastq PathToIMP3/test/Reads/test.mg.2.fastq"
  Metatranscriptomics: "PathToIMP3/test/Reads/test.mt.1.fastq PathToIMP3/test/Reads/test.mt.2.fastq"
sample: test
outputdir: "PathToIMP3/test/testrun_imp"
summarydir: "PathToIMP3/test/testrun_imp/summary"
summary_steps: "multiqc"
db_path: "PathToIMP3/databases"
pathofact:
  signalp: "PathToSignalP/signalp"
```
Repalce `PathToIMP3` and `PathToSignalP` with your paths.