# summary takes care of collecting rules and archiving data
localrules: SUMMARY, VIS

# include rules for the Report part of the workflow
sum_target = []
if "stats" in SUM_STEPS:
    include:
        '../../Summary/stats_gather.smk'
    sum_target.append("status/stats.done")
    if not "preprocessing" in IMP_STEPS and (MG or MT or MGaln or MTaln):
        include:
            '../../Summary/reads-count.smk'
        include:
            '../../Summary/stats_preprocessing.smk'
    if not "assembly" in IMP_STEPS and CONTIGS:
        include:
            '../../Assembly/common/mapping.smk'
        include:
            '../../Assembly/common/bwa.index.smk'
        include:
            '../../Assembly/common/contig.smk'
        include:
            '../../Assembly/common/contig-length.smk'
        include:
            '../../Assembly/common/contig-depth.smk'
        include:
            '../../Summary/stats_assembly.smk'

if "DB" in SUM_STEPS:
    include:
        '../../Summary/mongoDB.smk'
    sum_target.append("status/DB_filling.done")


if "vis" in SUM_STEPS:
    sum_target.append("status/visualization.done")
    vis_targets = []
    if "stats" in SUM_STEPS:
        include:
            '../../Summary/vis_stats.smk'
        vis_targets.append("status/stats.visualized")
    rule VIS:
        input:
            vis_targets
        output:
            touch('status/visualization.done')


# master command
rule SUMMARY:
    input:
        sum_target
    output:
        touch('status/summary.done')
