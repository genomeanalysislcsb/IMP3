# include global functions
include:
    "rules/function.definitions"

# include configuration file
include:
    "rules/ini/config"

# set working directory and dump output
workdir:
    OUTPUTDIR

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))
yaml.dump(config, open_output('full.config.yaml'), allow_unicode=True,default_flow_style=False)

# include rules for the workflow based on the input parameters
include:
    "rules/data.input.rules"

# include workflows for preprocessing, assembly, analysis, binning and taxonomy steps
# INTEGRATIVE MG-MT workflow
if MG and MT:
    if 'preprocessing' in IMP_STEPS:
        include:
            "workflows/integrative/Preprocessing"

    if 'assembly' in IMP_STEPS:
        if config['assembly']['hybrid']:
            include:
                "workflows/integrative/Assembly"
        else:
            include:
                "workflows/single_omics/mg/Assembly"

    if 'analysis' in IMP_STEPS:
        include:
            "workflows/integrative/Analysis"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflows/integrative/Taxonomy"

    if 'binning' in IMP_STEPS:
        include:
            "workflows/integrative/Binning"

    if 'summarize' in IMP_STEPS:
        include:
            "workflows/integrative/Summary"


# Single omics MG workflow
elif MG:
    if 'preprocessing' in IMP_STEPS:
        include:
            "workflows/single_omics/mg/Preprocessing"

    if 'assembly' in IMP_STEPS:
        include:
            "workflows/single_omics/mg/Assembly"

    if 'analysis' in IMP_STEPS:
        include:
            "workflows/single_omics/mg/Analysis"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflows/single_omics/mg/Taxonomy"

    if 'binning' in IMP_STEPS:
        include:
            "workflows/single_omics/mg/Binning"
# MT-only workflow
elif MT:
    if 'preprocessing' in IMP_STEPS:
        include:
            "workflows/single_omics/mt/Preprocessing"

    if 'assembly' in IMP_STEPS:
        include:
            "workflows/single_omics/mt/Assembly"

    if 'analysis' in IMP_STEPS:
        include:
            "workflows/single_omics/mt/Analysis"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflows/single_omics/mt/Taxonomy"

    if 'binning' in IMP_STEPS:
        include:
            "workflows/single_omics/mt/Binning"
#starting from aligned reads:
elif MGaln and MTaln:
    if 'analysis' in IMP_STEPS:
        include:
            "workflows/integrative/Analysis"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflows/integrative/Taxonomy"

    if 'binning' in IMP_STEPS:
        include:
            "workflows/integrative/Binning"

elif MGaln:
    if 'analysis' in IMP_STEPS:
        include:
            "workflows/single_omics/mg/Analysis"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflows/single_omics/mg/Taxonomy"

    if 'binning' in IMP_STEPS:
        include:
            "workflows/single_omics/mg/Binning"

elif MTaln:
    if 'analysis' in IMP_STEPS:
        include:
            "workflows/single_omics/mt/Analysis"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflows/single_omics/mt/Taxonomy"

    if 'binning' in IMP_STEPS:
        include:
            "workflows/single_omics/mt/Binning"
else:
    raise Exception('No input data.')

# master command summary and archiving
inputs_all = ['workflow.done']
if "summary" in IMP_STEPS:
    inputs_all.append('summary.done')
    include:
        "workflows/common/Summary"
#                 #if config['compress_level'] != "none" and config['compress_level'] != "":
#                 #    inputs_all.append('archiving.done')

localrules: workflow, ALL
# master command
rule ALL:
    input:
        inputs_all
    output:
        touch('all.done')

# master command workflow
rule workflow:
    input:
        workflow_ctrl
    output:
        touch('workflow.done')

# master command summary and archiving
inputs_all = ['workflow.done']
if "summary" in IMP_STEPS:
    inputs_all.append('summary.done')
    include:
        "workflows/common/Summary"
#if config['compress_level'] != "none" and config['compress_level'] != "":
#    inputs_all.append('archiving.done')

