include:
    '../Preprocessing/fastqc.smk'

include:
    '../Preprocessing/reads-count.smk'

if has_mg() and has_sr("mg", "pe"):
    include:
        '../Preprocessing/nonpareil.smk'

include:
    '../Preprocessing/archive.smk'

# PE
if (has_mg() or has_mt()) and (has_sr("mg", "pe") or has_sr("mt", "pe")):
    include:
        "../Preprocessing/trimming.smk"
    include:
        "../Preprocessing/rna.filtering.smk"
    include:
        "../Preprocessing/filtering.smk"

    rule PREPROCESSING_READY:
        input:
            expand(
                "Stats/{omic}/{omic}.{reads}_fastqc.{ext}",
                omic=get_omics(), reads=get_read_types_raw("mg"), ext=["zip", "html"]
            ),
            expand(
                "Stats/{omic}/{omic}.{reads}.preprocessed_fastqc.{ext}",
                omic=get_omics(), reads=get_read_types_proc("mg"), ext=["zip", "html"]
            ),
            expand(
                "Stats/{omic}/{omic}.read_counts.txt",
                omic=get_omics()
            ),
            expand(
                "Stats/mg/mg.{reads}_nonpareil.diversity.txt",
                reads=["r1", "r2"]
            ) if has_sr("mg", "pe") else [],
        output:
            touch("status/preprocessing.ready")
    
    rule PREPROCESSING:
        input:
            "status/preprocessing.ready",
            expand("{fname}.gz", fname=get_read_files_all(True)),
        output:
            touch("status/preprocessing.done")

# MG, PE
elif get_omics() == ("mg") and has_sr("mg", "pe"):
    # include:
    #     '../../../Preprocessing/reads-count.smk'
    # include:
    #     "../../../Preprocessing/trimming.smk"
    # include:
    #     "../../../Preprocessing/nonpareil.smk"
    # if PREPROCESSING_FILTERING:
    #     include:
    #         "../../../Preprocessing/filtering.smk"
    # else:
    #     include:
    #         "../../../Preprocessing/no-filtering.smk"
    # master command
    rule PREPROCESSING:
        input:
            # 'Preprocessing/mg.r1.preprocessed.fq',
            # 'Preprocessing/mg.r2.preprocessed.fq',
            # 'Preprocessing/mg.se.preprocessed.fq',
            # 'Stats/mg/mg.r1.preprocessed_fastqc.html',
            # 'Stats/mg/mg.r2.preprocessed_fastqc.html',
            # 'Stats/mg/mg.se.preprocessed_fastqc.html',
            # 'Stats/mg/mg.r1.preprocessed_fastqc.zip',
            # 'Stats/mg/mg.r2.preprocessed_fastqc.zip',
            # 'Stats/mg/mg.se.preprocessed_fastqc.zip',
            # 'Stats/mg/mg.read_counts.txt',
            # 'status/preprocessing_mg.reads.compressed',
            # expand('Stats/mg/mg.{read}_nonpareil.diversity.txt',read=["r1","r2"])
        output:
            touch('status/preprocessing.done')
elif get_omics() == ('mg') and not has_sr("mg", "pe"):
    # include:
    #     '../../../Preprocessing/reads-count_SE.smk'
    # include:
    #     "../../../Preprocessing/trimming_SE.smk"
    # if PREPROCESSING_FILTERING:
    #     include:
    #         "../../../Preprocessing/filtering_SE.smk"
    # else:
    #     include:
    #         "../../../Preprocessing/no-filtering_SE.smk"
    # master command
    rule PREPROCESSING:
        input:
            # expand([
            # 'Preprocessing/mg.r1.preprocessed.fq',
            # 'Preprocessing/mg.r2.preprocessed.fq',
            # 'Preprocessing/mg.se.preprocessed.fq',
            # "Stats/mg/mg.se.preprocessed_fastqc.zip",
            # "Stats/mg/mg.se.preprocessed_fastqc.html",
            # "Stats/mg/mg.read_counts.txt",
            # 'status/preprocessing_mg.reads.compressed'], type=TYPES)
        output:
            touch('status/preprocessing.done')
elif get_omics() == ("mt") and has_sr("mt", "pe"):
    # include:
    #     '../../../Preprocessing/reads-count.smk'
    # include:
    #     "../../../Preprocessing/trimming.smk"
    # include:
    #     "../../../Preprocessing/rna.filtering.smk"
    # if PREPROCESSING_FILTERING:
    #     include:
    #         "../../../Preprocessing/filtering.smk"
    # else:
    #     include:
    #         "../../../Preprocessing/no-filtering.smk"
    # master command
    rule PREPROCESSING:
        input:
            # 'Preprocessing/mt.r1.preprocessed.fq',
            # 'Preprocessing/mt.r2.preprocessed.fq',
            # 'Preprocessing/mt.se.preprocessed.fq',
            # 'Stats/mt/mt.r1.preprocessed_fastqc.html',
            # 'Stats/mt/mt.r2.preprocessed_fastqc.html',
            # 'Stats/mt/mt.se.preprocessed_fastqc.html',
            # 'Stats/mt/mt.r1.preprocessed_fastqc.zip',
            # 'Stats/mt/mt.r2.preprocessed_fastqc.zip',
            # 'Stats/mt/mt.se.preprocessed_fastqc.zip',
            # 'Stats/mt/mt.read_counts.txt',
            # 'status/preprocessing_mt.reads.compressed'
        output:
            touch('status/preprocessing.done')
