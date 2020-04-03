# include global functions
include:
    "workflow/rules/function.definitions.smk"

# include configuration file
include:
    "workflow/rules/ini/config.smk"

# set working directory and dump output
workdir:
    OUTPUTDIR

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))
yaml.dump(config, open_output('sample.config.yaml'), allow_unicode=True,default_flow_style=False)

Path("status").mkdir(exist_ok=True)

# include rules for the workflow based on the input parameters
include:
    "workflow/rules/data.input.smk"

# include workflows for preprocessing, assembly, analysis, binning and taxonomy steps
# INTEGRATIVE MG-MT workflow
if MG and MT:
    if 'preprocessing' in IMP_STEPS:
        include:
            "workflow/rules/modules/integrative/Preprocessing.smk"

    if 'assembly' in IMP_STEPS:
        if config['assembly']['hybrid']:
            include:
                "workflow/rules/modules/integrative/Assembly.smk"
        else:
            include:
                "workflow/rules/modules/single_omics/mg/Assembly.smk"

    if 'analysis' in IMP_STEPS:
        include:
            "workflow/rules/modules/integrative/Analysis.smk"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflow/rules/modules/integrative/Taxonomy.smk"

    if 'binning' in IMP_STEPS:
        include:
            "workflow/rules/modules/integrative/Binning.smk"

    if 'summary' in IMP_STEPS:
        include:
            "workflow/rules/modules/common/Summary.smk"


# Single omics MG workflow
elif MG:
    if 'preprocessing' in IMP_STEPS:
        if len(MG) == 2:
            include:
                "workflow/rules/modules/single_omics/mg/Preprocessing.smk"
        if len(MG) == 1:
            include:
                "workflow/rules/modules/single_omics/mg/PreprocessingSE.smk"

    if 'assembly' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mg/Assembly.smk"

    if 'analysis' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mg/Analysis.smk"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mg/Taxonomy.smk"

    if 'binning' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mg/Binning.smk"
# MT-only workflow
elif MT:
    if 'preprocessing' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Preprocessing.smk"

    if 'assembly' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Assembly.smk"

    if 'analysis' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Analysis.smk"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Taxonomy.smk"

    if 'binning' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Binning.smk"
#starting from aligned reads:
elif MGaln and MTaln:
    if 'analysis' in IMP_STEPS:
        include:
            "workflow/rules/modules/integrative/Analysis.smk"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflow/rules/modules/integrative/Taxonomy.smk"

    if 'binning' in IMP_STEPS:
        include:
            "workflow/rules/modules/integrative/Binning.smk"

elif MGaln:
    if 'analysis' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mg/Analysis.smk"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mg/Taxonomy.smk"

    if 'binning' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mg/Binning.smk"

elif MTaln:
    if 'analysis' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Analysis.smk"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Taxonomy.smk"

    if 'binning' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Binning.smk"
else:
    raise Exception('No input data.')

# master command summary and archiving
inputs_all = ['status/inputs.done','status/workflow.done']
if "summary" in IMP_STEPS:
    inputs_all.append('status/summary.done')
    include:
        "workflow/rules/modules/common/Summary.smk"


# master command
rule ALL:
    input:
        inputs_all
    output:
        touch('status/all.done')


localrules: workflow, ALL
# master command workflow
rule workflow:
    input:
        workflow_ctrl
    output:
        touch('status/workflow.done')

