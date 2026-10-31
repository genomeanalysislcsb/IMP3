# include rules for the Preprocessing part of the workflow

include:
    "../../../Preprocessing/fastqc.smk"

include:
    '../../../Preprocessing/reads-count.smk'

include:
    "../../../Preprocessing/nonpareil.smk"

# trimming rules
include:
    "../../../Preprocessing/trimming.smk"

# include filtering rule or not
if PREPROCESSING_FILTERING:
    include:
        "../../../Preprocessing/filtering.smk"
else:
    include:
        "../../../Preprocessing/no-filtering.smk"

# master command
rule PREPROCESSING:
    input:
        'Preprocessing/mg.r1.preprocessed.fq.gz',
        'Preprocessing/mg.r2.preprocessed.fq.gz',
        'Preprocessing/mg.se.preprocessed.fq.gz',
        'Stats/mg/mg.r1.preprocessed_fastqc.html',
        'Stats/mg/mg.r2.preprocessed_fastqc.html',
        'Stats/mg/mg.se.preprocessed_fastqc.html',
        'Stats/mg/mg.r1.preprocessed_fastqc.zip',
        'Stats/mg/mg.r2.preprocessed_fastqc.zip',
        'Stats/mg/mg.se.preprocessed_fastqc.zip',
        'Stats/mg/mg.read_stats.tsv',
        # 'status/preprocessing_mg.reads.compressed',
        expand('Stats/mg/mg.{read}_nonpareil.diversity.txt',read=["r1","r2"])
    output:
        touch('status/preprocessing.done')
