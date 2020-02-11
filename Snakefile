# include configuration file
include:
    "rules/ini/config"


# define the data types used and the assembly
if MG and MT and config['assembly']['hybrid']:
    TYPES = ['mg', 'mt']
    ASS = 'mgmt'
elif MG and MT:
    TYPES = ['mg', 'mt']
    ASS = 'mg'
elif MG:
    TYPES = ['mg']
    ASS = 'mg'
elif MT:
    TYPES = ['mt']
    ASS = 'mt'
elif MGaln and MTaln:
    TYPES = ['mg', 'mt']
    ASS = 'mg'
elif MGaln:
    TYPES = ['mg']
    ASS = 'mg'
elif MTaln:
    TYPES = ['mt']
    ASS = 'mt'
elif CONTIGS and MG and MT:
    TYPES = ['mg', 'mt']
    ASS = 'mg'
elif CONTIGS and MG:
    TYPES = ['mg']
    ASS = 'mg'
elif CONTIGS and MT:
    TYPES = ['mt']
    ASS = 'mg'

workdir:
    OUTPUTDIR

f = open(OUTPUTDIR+'/full.config.yaml', 'w+')
yaml.dump(config, f, allow_unicode=True,default_flow_style=False)


# include rules for the workflow based on the input parameters
include:
    "rules/data.input.rules"

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

inputs = []
if 'preprocessing' in IMP_STEPS:
    inputs.append('preprocessing.done')
if 'assembly' in IMP_STEPS:
    inputs.append('assembly.done')
if 'analysis' in IMP_STEPS:
    inputs.append('analysis.done')
if 'taxonomy' in IMP_STEPS:
    inputs.append('taxonomy.done')
if 'binning' in IMP_STEPS:
    inputs.append('binning.done')

# master command
rule ALL:
    input:
        inputs
    output:
        touch('workflow.done')
