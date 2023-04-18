# include rules for the Preprocessing part of the workflow

include:
    "../../../Preprocessing/fastqc.smk"
include:
    '../../../Preprocessing/reads-count.smk'
include:
    "../../../Preprocessing/trimming.smk"
include:
    "../../../Preprocessing/rna.filtering.smk"

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
        'Preprocessing/mt.r1.preprocessed.fq',
        'Preprocessing/mt.r2.preprocessed.fq',
        'Preprocessing/mt.se.preprocessed.fq',
        'Stats/mt/mt.r1.preprocessed_fastqc.html',
        'Stats/mt/mt.r2.preprocessed_fastqc.html',
        'Stats/mt/mt.se.preprocessed_fastqc.html',
        'Stats/mt/mt.r1.preprocessed_fastqc.zip',
        'Stats/mt/mt.r2.preprocessed_fastqc.zip',
        'Stats/mt/mt.se.preprocessed_fastqc.zip',
        'Stats/mt/mt.read_counts.txt',
        # 'status/preprocessing_mt.reads.compressed'
    output:
        touch('status/preprocessing.done')
