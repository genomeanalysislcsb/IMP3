binctrl = []
if config['dereplication']['simple_derep']:
    include:
        "../../Multi/Dereplicate/dRep.simple.smk"
    include:
            "../../Multi/Dereplicate/dRep_sum.simple.smk"
    binctrl.append("Dereplication/dereplicated_bin_stats.tsv")
if config['dereplication']['cross_mapping_rebinning']['do']:
    include:
        "../../Multi/Dereplicate/cross_mapping.smk"
    include:
        "../../Multi/Dereplicate/MetaBAT.smk"
    include:
        "../../Multi/Dereplicate/MaxBin.smk"
    include:
        "../../Multi/Dereplicate/dastool.smk"
    include:
        "../../Multi/Dereplicate/grid.smk"
    include:
        "../../Multi/Dereplicate/dRep.multi.smk"
    binctrl.append("Dereplication/rebinned_bin_stats.tsv")
    if config['dereplication']['GTDBtk']:
        include:
            "../../Multi/Dereplicate/gtdbtk_multi.smk"
    else:
        include:
            "../../Multi/Dereplicate/dRep_sum.multi.smk"


# master command
rule all_derep:
    input:
        binctrl
    output:
        touch('status/dereplicate.done')


