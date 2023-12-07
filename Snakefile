import pandas as pd
import numpy as np

# include global functions
ROOTDIR = srcdir(".")

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
yaml.dump(config, open_output('sample.config.yaml'), allow_unicode=True, default_flow_style=False)

Path("status").mkdir(exist_ok=True)

sample_input_tsv = config['sample_input_tsv']

# Samples
SAMPLES = pd.read_csv(sample_input_tsv, sep="\t", header=0, index_col="SampleID").replace({pd.NA: None, 'NA': None, np.nan: None})

# List of all sample IDs
sample_ids = list(SAMPLES.index)

# List of all file types
file_types = list(SAMPLES.columns)

# List of all file names
input_read_file_names = ["mg.r1.fq", "mg.r2.fq", "mg.se.fq", "mt.r1.fq", "mt.r2.fq", "mt.se.fq", "assembly.fna"]

# Dict linking file types to input file names
# (we use this in get_ori_files() to get the correct file name for each file type)
input_read_file_names_dict = {name: type for name, type in zip(input_read_file_names, file_types)}

# Reverse dict to get file type from file name
input_read_file_names_dict_rev = {type: name for name, type in input_read_file_names_dict.items()}

### TEST FUNCTIONS ###

def get_read_out_paths_and_wc(filt_samples, filt_types):
    filt_samples, filtered_file_reads, filtered_file_omes = zip(*[[sample] + file_type.split(".")[:-1] for sample, file_type in zip(filt_samples, filt_types) if file_type.endswith(".fq")])

    # Build expected file copy output paths
    # (notice the `zip` to not create all combinations of sample IDs and file types)
    valid_file_paths = expand(os.path.join(OUTPUTDIR, "{sid}", "Preprocessing", "{ome}.{read}.fq"), zip,
                              sid=filt_samples, ome=filtered_file_reads, read=filtered_file_omes)

    checked_file_paths = expand(os.path.join(OUTPUTDIR, "{sid}", "status", "{ome}.{read}.fq.checked"), zip,
                            sid=filt_samples, ome=filtered_file_reads, read=filtered_file_omes)
    return valid_file_paths, checked_file_paths

# TODO: Finish assembly output paths function
def get_assembly_out_paths_and_wc(filt_samples, filt_types):
    # Build expected file copy output paths, matching the sample ID and exisitng assembly file name
    filt_samples, filtered_file_assemblies = zip(*[(sample, file_type.split(".")[0]) for sample, file_type in zip(filt_samples, filt_types) if file_type == "assembly.fna"])

    valid_file_paths = expand(os.path.join(OUTPUTDIR, "{sid}", "Assembly", "{ass}.fna"), zip,
                              sid=filt_samples, ass=filtered_file_assemblies)

    checked_file_paths = expand(os.path.join(OUTPUTDIR, "{sid}", "status", "{ass}.fna.checked"), zip,
                            sid=filt_samples, ass=filtered_file_assemblies)

    return valid_file_paths, checked_file_paths

def get_files_to_copy(samples_df):
    # Stack the DataFrame and reset index to get easy access to rows and columns
    stacked_df = samples_df.stack().reset_index()

    # Extract sample(rows) and file types columns, excluding rows with 'None' values
    filtered_samples, filtered_file_types = zip(*[(row, input_read_file_names_dict_rev.get(col)) 
                                                   for row, col, val in stacked_df.values if val])

    valid_read_file_paths, checked_read_file_paths = get_read_out_paths_and_wc(filtered_samples, filtered_file_types)
    valid_assembly_file_paths, checked_assembly_file_paths = get_assembly_out_paths_and_wc(filtered_samples, filtered_file_types)

    # print('filtered_samples:', filtered_samples)
    # print('valid_read_file_paths:', valid_read_file_paths)
    # print('checked_read_file_paths:', checked_read_file_paths)
    # print('valid_assembly_file_paths', valid_assembly_file_paths)
    # print('checked_assembly_file_paths', checked_assembly_file_paths)

    combined_checked_paths = checked_read_file_paths + checked_assembly_file_paths

    print('combined_checked_paths', combined_checked_paths)
    print(all([isinstance(path, str) for path in combined_checked_paths]))

    return combined_checked_paths

def get_ori_files(wildcards):
    # Get sample ID based on what is defined in rule All based on whats in the sample_input_tsv (SAMPLES)
    sample_id = wildcards.sid

    # Get file type in the same manner
    file_type = input_read_file_names_dict.get(wildcards.file)

    # Get file path based on sample ID and file type
    sample_id_file = SAMPLES.loc[sample_id, file_type]

    return sample_id_file
    
###########################

### RULES #################
localrules: ALL
####
rule ALL:
    input:
        get_files_to_copy(SAMPLES)
####
rule prepare_input_data_copy:
    input:
        get_ori_files # --> ['/path/to/mg1.fq', '/path/to/mg2.fq', ...]
    output:
        os.path.join(OUTPUTDIR, "{sid}", "Preprocessing", "{file}")
    threads:
        getThreads(5)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    wildcard_constraints:
        sid = "|".join(sample_ids)
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "Prepare input: decompress/copy"
    script:
        os.path.join(SRCDIR, "prepare_input_files.py")

rule check_read_r_input:
    input:
        expand(os.path.join(OUTPUTDIR, "{{sid}}", "Preprocessing", "{{ome}}.{read}.fq"), read=["r1","r2"])
    output:
        expand(os.path.join(OUTPUTDIR, "{{sid}}", "status", "{{ome}}.{read}.fq.checked"), read=["r1","r2"])
    priority:
        10
    message:
        "Checking {wildcards.ome} read input."
    shell:
        """
        if [ "$(wc -l < {input[0]})" -eq "$(wc -l < {input[1]})" ]; then
          touch {output}
        else
          echo "ERROR: Paired read files don't have the same number of reads." && exit 1
        fi
        """

rule check_read_se_input:
    input:
        expand(os.path.join(OUTPUTDIR, "{{sid}}", "Preprocessing", "{{ome}}.{read}.fq"), read=["se"])
    output:
        expand(os.path.join(OUTPUTDIR, "{{sid}}", "status", "{{ome}}.{read}.fq.checked"), read=["se"])
    message:
        "Checking single end {wildcards.ome} read input."
    priority:
        10
    shell:
        """
        if [ -s {input} ]; then
          touch {output}
        else
          echo "ERROR: {input} is empty." && exit 1
        fi
        """

# TODO: Add check for empty assembly files
rule check_assembly_input:
    input:
        "Assembly/{ass}.assembly.merged.fa"
    output:
        "status/{ass}.assembly.merged.fa.checked"
    priority:
        10
    message:
        "Checking assembly."
    shell:
        """
        if [ -s {input} ]; then
          touch {output}
        else
          echo "ERROR: Assembly is empty." && exit 1
        fi
        """

###########################
