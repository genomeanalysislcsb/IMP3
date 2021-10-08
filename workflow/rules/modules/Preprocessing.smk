if has_sr("mg", "pe") or has_sr("mt", "pe"):
    include:
        "../Preprocessing/trim.smk"
    include:
        "../Preprocessing/filter.smk"
else:
    include:
        "../Preprocessing/trim_se.smk"
    include:
        "../Preprocessing/filter_se.smk"

if has_mt():
    include:
        "../Preprocessing/filter_rrna.smk"

include:
    '../Preprocessing/fastqc.smk'

if has_mg() and has_sr("mg", "pe"):
    include:
        '../Preprocessing/nonpareil.smk'

include:
    '../Preprocessing/misc.smk'

rule PREPROCESSING_READY:
    input:
        expand(
            "Stats/{omic}/{omic}.{reads}_fastqc.{ext}",
            omic=get_omics(), reads=get_read_types_raw(), ext=["zip", "html"]
        ),
        expand(
            "Stats/{omic}/{omic}.{reads}.preprocessed_fastqc.{ext}",
            omic=get_omics(), reads=get_read_types_proc(), ext=["zip", "html"]
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
