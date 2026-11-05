# include global functions
ROOTDIR = srcdir('')

# print(f'workflow.cores, workflow.nodes: {workflow.cores}, {workflow.nodes}')
# print(f'workflow.conda_prefix: {workflow.conda_prefix}, {os.environ["CONDA_PREFIX"]}')

include:
    "workflow/rules/function.definitions.smk"

# include configuration file
include:
    "workflow/rules/ini/config.smk"

# Optional: run every rule in the IMP3 container (Apptainer/Singularity image with all conda environments, see
# docs/portable.md). Set `container: /path/to/imp3.sif` in the config and run with --use-conda --use-singularity.
if config.get('container'):
    containerized:
        config['container']

# set working directory and dump output
workdir:
    OUTPUTDIR

# If we have `status/all.done` exit without error and with success message
if os.path.exists('status/all.done'):
    import sys
    sys.exit('IMP3 already finished successfully.')

# Break if we have `status/workflow.done` + `status/stats.done`, since something wnt wrong but everything is already done
required_files = ['status/workflow.done', 'status/stats.done', 'status/compress_text_files.done', 'status/compress_dirs.done', 'status/summary.done']
if all(os.path.exists(file) for file in required_files) and not os.path.exists('status/all.done'):
    # Make all.done file
    Path('status/all.done').touch()
    raise Exception('Workflow already finished, but final status file is missing.')


yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))
yaml.dump(config, open_output('sample.config.yaml'), allow_unicode=True, default_flow_style=False)

# Ensure some directories exist
Path("status").mkdir(exist_ok=True)
Path("status/logs").mkdir(exist_ok=True)
Path("logs_cluster_jobs").mkdir(exist_ok=True)

# include rules for the workflow based on the input parameters
include:
    "workflow/rules/data.input.smk"

# TODO: Add user param to specify interleaved read input and add a rule to check for it and deiniterleave it if necessary

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

# Single omics MG workflow
elif MG:
    if 'preprocessing' in IMP_STEPS:
        if len(MG) == 2:
            include:
                "workflow/rules/modules/single_omics/mg/Preprocessing.smk"
        if len(MG) == 1:
            include:
                "workflow/rules/modules/single_omics/mg/Preprocessing_SE.smk"

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
        if len(MT) == 2:
            include:
                "workflow/rules/modules/single_omics/mt/Preprocessing.smk"
        if len(MT) == 1:
            include:
                "workflow/rules/modules/single_omics/mt/Preprocessing_SE.smk"

    if 'assembly' in IMP_STEPS:
        if len(MT) == 2:
            include:
                "workflow/rules/modules/single_omics/mt/Assembly.smk"
        if len(MT) == 1:
            include:
                "workflow/rules/modules/single_omics/mt/Assembly_SE.smk"
    
    if 'analysis' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Analysis.smk"

    if 'taxonomy' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Taxonomy.smk"

    if 'binning' in IMP_STEPS:
        include:
            "workflow/rules/modules/single_omics/mt/Binning.smk"
# Starting from aligned reads:
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

# Clean up at the end
if EMAIL == "":
    onsuccess:
        shell("""
              tar -czf .snakemake.tar.gz .snakemake &>> status/logs/cleanup.log
              rm -rf {config[tmp_dir]} &>> status/logs/cleanup.log
              rm -rf .snakemake &>> status/logs/cleanup.log
              """)
else:
    onsuccess:
        shell("""
              echo "$(date) {config[sessionName]}" | mail -s "IMP3 finished" {EMAIL}
              tar -czf .snakemake.tar.gz .snakemake &>> status/logs/cleanup.log
              rm -rf {config[tmp_dir]} &>> status/logs/cleanup.log
              rm -rf .snakemake &>> status/logs/cleanup.log
              """)
    onerror:
        shell('echo "$(date) {config[sessionName]}" | mail -s "IMP3 exited with error" {EMAIL}')
    onstart:
        shell('echo "$(date) {config[sessionName]}" | mail -s "IMP3 started" {EMAIL}')

# Input for archiving
text_files_to_compress_analysis = []
if 'analysis' in IMP_STEPS:
    text_files_to_compress_analysis.append("Analysis/annotation/annotation_CDS_RNA_hmms.gff")
    # From featureCount.smk

    text_files_to_compress_analysis.append(aggregate_mantis_fc if config["annotation"] == "mantis" else [])

text_files_to_compress_binning = []
if 'binning' in IMP_STEPS:
    text_files_to_compress_binning.append("Binning/MetaWrap/bins.stats")

# Master command summary and archiving
inputs_all = ['status/inputs.done', 'status/workflow.done', 'status/compress_dirs.done', 'status/compress_text_files.done']

if 'preprocessing' in IMP_STEPS:
    inputs_all.append(expand('status/preprocessing_{type}.reads.compressed', type=TYPES))

if 'assembly' in IMP_STEPS:
    inputs_all.append("Assembly/intermediary.tar.gz")

if 'analysis' in IMP_STEPS:
    # Add the gzipped versions of the text files to inputs_all (only if its not a function (e.g. checkpint function))
    inputs_all.extend([f"{file}.gz" for file in text_files_to_compress_analysis if file != aggregate_mantis_fc])
    # We dont actually need to add this to the inputs_all list, as it is already included in the text_files_to_compress_analysis list and jsut causes trouble
    # inputs_all.append(aggregate_mantis_fc_gz if config["annotation"] == "mantis" else [])

if 'binning' in IMP_STEPS:
    inputs_all.extend([f"{file}.gz" for file in text_files_to_compress_binning])

if "summary" in IMP_STEPS:
    inputs_all.append('status/summary.done')
    include:
        "workflow/rules/modules/common/Summary.smk"


# if config['proteomics']['add_kraken_taxa_proteomes']:
#     inputs_all.append('status/tarball_taxa_proteins.done')

# List of text files to compress
text_files_to_compress = text_files_to_compress_analysis + text_files_to_compress_binning

localrules: workflow, ALL

# Master command
rule ALL:
    input:
        inputs_all
    output:
        touch('status/all.done')

# Rule for the main workflow
rule workflow:
    input:
        workflow_ctrl
    output:
        touch('status/workflow.done')

# Rule to compress directories and remove originals
rule compress_dirs:
    input:
        'status/workflow.done',
        'Stats/all_stats.tsv' if 'summary' in IMP_STEPS else [],
        'Stats/all_stats.Rdata' if 'summary' in IMP_STEPS else [],
        'Stats/all_bin_stats.tsv' if 'binning' in IMP_STEPS else [],
        'Assembly/intermediary.tar.gz' if 'assembly' in IMP_STEPS else [],
        expand('Preprocessing/{type}.{read}.preprocessed.fq.gz', type=TYPES, read=['r1', 'r2', 'se']) if 'preprocessing' in IMP_STEPS else []
    output:
        touch('status/compress_dirs.done')
    params:
        logs_dir = 'logs',
        cluster_logs_dir = 'logs_cluster_jobs'
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        4
    conda:
        os.path.join(ENVDIR, "IMP_compress.yaml")
    log:
        'status/logs/compress_dirs.log'
    message:
        "Compressing directories: logs, logs_cluster_jobs"
    shell:
        """
        if [ -d '{params.logs_dir}' ]; then
            echo "Compressing {params.logs_dir}" >> {log} 2>&1
            tar -czf {params.logs_dir}.tar.gz {params.logs_dir} >> {log} 2>&1
            rm -rf {params.logs_dir} >> {log} 2>&1
        fi
        if [ -d '{params.cluster_logs_dir}' ]; then
            echo "Compressing {params.cluster_logs_dir}" >> {log} 2>&1
            tar -czf {params.cluster_logs_dir}.tar.gz {params.cluster_logs_dir} >> {log} 2>&1
            rm -rf {params.cluster_logs_dir} >> {log} 2>&1
        fi
        touch {output} >> {log} 2>&1
        """

# Rule to compress text-based files larger than 1KB
rule compress_text_files:
    input:
        'status/compress_dirs.done',
        text_files_to_compress,
    output:
        'Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz',
        'Binning/MetaWrap/bins.stats.gz',
        touch('status/compress_text_files.done')
    params:
        min_size = 0  # in bytes
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(16)
    conda:
        os.path.join(ENVDIR, "IMP_compress.yaml")
    log:
        'status/logs/compress_text_files.log'
    message:
        "Compressing text files larger than {params.min_size} bytes"
    shell:
        """
        echo "Compressing text files larger than {params.min_size} bytes ..." >> {log} 2>&1
        MIN_SIZE={params.min_size}
        find . -type f ! -name "*.gz" ! -name "*.tar.gz" ! -path "./.snakemake/*" -size +${{MIN_SIZE}}c -exec file {{}} + | grep -E 'text|ASCII' | cut -d: -f1 | parallel -j {threads} gzip -f >> {log} 2>&1
        touch {output} >> {log} 2>&1
        echo "Done." >> {log} 2>&1
        """

# if config['proteomics']['add_kraken_taxa_proteomes']:
#     rule tarball_taxa_proteins:
#         input:
#             'status/compress_text_files.done'
#         output:
#             check=touch('status/tarball_taxa_proteins.done')
#         params:
#             tax_prot_dl_dir='Analysis/annotation/proteomics/taxa_proteins/download',
#             target='files',
#             tarball='Analysis/annotation/proteomics/taxa_proteins/download/files.tar',
#         resources:
#             runtime = VERY_SHORT_RUNTIME,
#             mem_mb = MAIN_MEMCORE_MB
#         threads:
#             getThreads(4)
#         conda:
#             os.path.join(ENVDIR, "IMP_compress.yaml")
#         log:
#             'status/logs/tarball_taxa_proteins.log'
#         message:
#             "Creating tar archive of reference taxa proteome files and deleting the original directory"
#         shell:
#             """
#             echo "Tarballing {params.target} ...">> {log}
#             tar -cf {params.tarball} -C {params.tax_prot_dl_dir} {params.target} >> {log} 2>&1
#             rm -rf {params.target} >> {log} 2>&1
#             touch {output.check} >> {log} 2>&1
#             echo "Done.">> {log}
#             """
