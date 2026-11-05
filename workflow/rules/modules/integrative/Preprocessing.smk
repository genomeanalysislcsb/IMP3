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
        'Preprocessing/{type}.r1.preprocessed.fq.gz',
        'Preprocessing/{type}.r2.preprocessed.fq.gz',
        'Preprocessing/{type}.se.preprocessed.fq.gz',
        "Stats/{type}/{type}.r1.preprocessed_fastqc.zip",
        "Stats/{type}/{type}.r2.preprocessed_fastqc.zip",
        "Stats/{type}/{type}.se.preprocessed_fastqc.zip",
        "Stats/{type}/{type}.r1.preprocessed_fastqc.html",
        "Stats/{type}/{type}.r2.preprocessed_fastqc.html",
        "Stats/{type}/{type}.se.preprocessed_fastqc.html",
        "Stats/{type}/{type}.read_stats.tsv"], type=TYPES),
        expand('Stats/mg/mg.{read}_nonpareil.diversity.txt',read=["r1","r2"])]

# metaT reads for the assembly (coverage normalization with BBNorm was removed)
MT_R1_ASSEMBLY = 'Preprocessing/mt.r1.preprocessed.fq.gz'
MT_R2_ASSEMBLY = 'Preprocessing/mt.r2.preprocessed.fq.gz'
MT_SE_ASSEMBLY = 'Preprocessing/mt.se.preprocessed.fq.gz'

# master command
rule PREPROCESSING:
    input:
        PREP_IN
    output:
        touch('status/preprocessing.done')
