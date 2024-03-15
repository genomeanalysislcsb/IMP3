# include rules for the Preprocessing part of the workflow

include:
    '../../Preprocessing/fastqc.smk'

include:
    '../../Preprocessing/reads-count.smk'

include:
    '../../Preprocessing/nonpareil.smk'

# trimming rules
include:
    "../../Preprocessing/trimming.smk"

# include filtering rule or not
if PREPROCESSING_FILTERING:
    include:
        "../../Preprocessing/filtering.smk"
else:
    include:
        "../../Preprocessing/no-filtering.smk"

# filtering rna for mt data
include:
    "../../Preprocessing/rna.filtering.smk"

# master rule input
PREP_IN = [expand([
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
        "Stats/{type}/{type}.r1.preprocessed_fastqc.zip",
        "Stats/{type}/{type}.r2.preprocessed_fastqc.zip",
        "Stats/{type}/{type}.se.preprocessed_fastqc.zip",
        "Stats/{type}/{type}.r1.preprocessed_fastqc.html",
        "Stats/{type}/{type}.r2.preprocessed_fastqc.html",
        "Stats/{type}/{type}.se.preprocessed_fastqc.html",
        "Stats/{type}/{type}.read_counts.txt"], type=TYPES),
        expand('Stats/mg/mg.{read}_nonpareil.diversity.txt',read=["r1","r2"])]

# include MT read coverafe normalization or not
if NORMALIZE_MT_READ_COVERAGE:
    include:
        "../../Preprocessing/normalize_coverage.smk"

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
