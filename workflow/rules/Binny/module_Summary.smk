# summary takes care of collecting rules and archiving data
localrules: SUMMARY, VIS

# include rules for the Report part of the workflow
sum_target = []
include:
    'stats_gather.smk'
sum_target.append("status/stats.done")
include:
    'contig.smk'
include:
    '../../Assembly/common/contig-length.smk'
include:
    'contig-depth.smk'

sum_target.append("status/visualization.done")
vis_targets = []
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
