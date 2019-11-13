# include configuration file
include:
    "rules/ini/config"

# define the data types used and the assembly
if MG and MT:
    TYPES = ['mg', 'mt']
    ASS = 'mgmt'
elif MG:
    TYPES = ['mg']
    ASS = 'mg'
elif MT:
    TYPES = ['mt']
    ASS = 'mt'


workdir:
    OUTPUTDIR

# include rules for the workflow based on the input parameters
include:
    "rules/data.input.rules"

# INTEGRATIVE MG-MT workflow
if MG and MT:
    if 'preprocessing' in IMP_STEPS:
        include:
            "workflows/integrative/Preprocessing"

    if 'assembly' in IMP_STEPS:
        include:
            "workflows/integrative/Assembly"

    if 'analysis' in IMP_STEPS:
        include:
            "workflows/integrative/Analysis"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflows/integrative/Taxonomy"

    if 'binning' in IMP_STEPS:
        include:
            "workflows/integrative/Binning"

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
