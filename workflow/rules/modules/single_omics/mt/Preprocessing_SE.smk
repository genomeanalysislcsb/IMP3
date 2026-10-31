# include rules for the Preprocessing part of the workflow

include:
    '../../../Preprocessing/fastqc_SE.smk'
include:
    '../../../Preprocessing/reads-count_SE.smk'
include:
    "../../../Preprocessing/trimming_SE.smk"
include:
    "../../../Preprocessing/rna.filtering_SE.smk"

# include filtering rule or not
if PREPROCESSING_FILTERING:
    include:
        "../../../Preprocessing/filtering_SE.smk"
else:
    include:
        "../../../Preprocessing/no-filtering_SE.smk"

# master rule input
PREP_IN = ["Stats/mt/mt.se_fastqc.html",
           "Stats/mt/mt.se_fastqc.zip",
           "Preprocessing/mt.se.preprocessed.fq.gz",
           "Stats/mt/mt.se.preprocessed_fastqc.zip",
           "Stats/mt/mt.se.preprocessed_fastqc.html",
           "Stats/mt/mt.read_stats.tsv",]


# include MT read coverafe normalization or not
if NORMALIZE_MT_READ_COVERAGE:
    include:
        "../../../Preprocessing/normalize_coverage.smk"

    MT_SE_ASSEMBLY = 'Preprocessing/mt.se.preprocessed.normalized.fq'

    PREP_IN.append(temp('Preprocessing/mt.se.preprocessed.normalized.fq'))
else:
    MT_SE_ASSEMBLY = 'Preprocessing/mt.se.preprocessed.fq.gz'

# master command
rule PREPROCESSING:
    input:
        PREP_IN
    output:
        touch('status/preprocessing.done')
