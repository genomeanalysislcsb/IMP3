# include rules for the Preprocessing part of the workflow

include:
    '../../../Preprocessing/fastqc_SE.smk'

include:
    '../../../Preprocessing/reads-count_SE.smk'

# trimming rules
include:
    "../../../Preprocessing/trimming_SE.smk"

# include filtering rule or not
if PREPROCESSING_FILTERING:
    include:
        "../../../Preprocessing/filtering_SE.smk"
else:
    include:
        "../../../Preprocessing/no-filtering_SE.smk"

# master command
rule PREPROCESSING:
    input:
        expand([
        'Preprocessing/mg.r1.preprocessed.fq',
        'Preprocessing/mg.r2.preprocessed.fq',
        'Preprocessing/mg.se.preprocessed.fq',
        "Stats/mg/mg.se.preprocessed_fastqc.zip",
        "Stats/mg/mg.se.preprocessed_fastqc.html",
        "Stats/mg/mg.read_stats.tsv",], type=TYPES)
        # 'status/preprocessing_mg.reads.compressed'], type=TYPES)
    output:
        touch('status/preprocessing.done')
