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
# filtering rna for mt data
include:
    "../../Preprocessing/rna.filtering.smk"
# include filtering rule or not
if PREPROCESSING_FILTERING:
    include:
        "../../Preprocessing/filtering.smk"
else:
    include:
        "../../Preprocessing/no-filtering.smk"

# master command
rule PREPROCESSING:
    input:
        expand([
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
        "Stats/{type}/{type}.r1.preprocessed_fastqc.zip",
        "Stats/{type}/{type}.r2.preprocessed_fastqc.zip",
        "Stats/{type}/{type}.se.preprocessed_fastqc.zip",
        "Stats/{type}/{type}.r1.preprocessed_fastqc.html",
        "Stats/{type}/{type}.r2.preprocessed_fastqc.html",
        "Stats/{type}/{type}.se.preprocessed_fastqc.html",
        "Stats/{type}/{type}.read_counts.txt",
        'status/preprocessing_{type}.reads.compressed'], type=TYPES),
        expand('Stats/mg/mg.{read}_nonpareil.diversity.txt',read=["r1","r2"])
    output:
        touch('status/preprocessing.done')
