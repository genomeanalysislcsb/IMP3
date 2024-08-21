# include rules for the Preprocessing part of the workflow

include:
    "../../../Preprocessing/fastqc_SE.smk"
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


# master rule input
PREP_IN = ['Preprocessing/mt.r1.preprocessed.fq',
        'Preprocessing/mt.r2.preprocessed.fq',
        'Preprocessing/mt.se.preprocessed.fq',
        'Stats/mt/mt.r1.preprocessed_fastqc.html',
        'Stats/mt/mt.r2.preprocessed_fastqc.html',
        'Stats/mt/mt.se.preprocessed_fastqc.html',
        'Stats/mt/mt.r1.preprocessed_fastqc.zip',
        'Stats/mt/mt.r2.preprocessed_fastqc.zip',
        'Stats/mt/mt.se.preprocessed_fastqc.zip',
        'Stats/mt/mt.read_stats.tsv']


# include MT read coverafe normalization or not
if NORMALIZE_MT_READ_COVERAGE:
    include:
        "../../../Preprocessing/normalize_coverage.smk"

    MT_R1_ASSEMBLY = 'Preprocessing/mt.r1.preprocessed.normalized.fq'
    MT_R2_ASSEMBLY = 'Preprocessing/mt.r2.preprocessed.normalized.fq'
    MT_SE_ASSEMBLY = 'Preprocessing/mt.se.preprocessed.normalized.fq'

    PREP_IN.append(temp('Preprocessing/mt.r1.preprocessed.normalized.fq'))
    PREP_IN.append(temp('Preprocessing/mt.r2.preprocessed.normalized.fq'))
    PREP_IN.append(temp('Preprocessing/mt.se.preprocessed.normalized.fq'))
else:
    MT_R1_ASSEMBLY = 'Preprocessing/mt.r1.preprocessed.fq'
    MT_R2_ASSEMBLY = 'Preprocessing/mt.r2.preprocessed.fq'
    MT_SE_ASSEMBLY = 'Preprocessing/mt.se.preprocessed.fq'

# master command
rule PREPROCESSING:
    input:
        PREP_IN
    output:
        touch('status/preprocessing.done')
