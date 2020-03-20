# include global functions
include:
    "workflow/rules/function.definitions.smk"

# include configuration file
include:
    "workflow/rules/ini/multi_config.smk"

# set working directory and dump output
workdir:
    OUTPUTDIR

yaml.dump(config, open_output('multi.config.yaml'), allow_unicode=True,default_flow_style=False)

# include workflows for data collection, catalogue, dRep and DB
if 'collate' in MULTI_STEPS:
    include:
        "workflow/rules/modules/multi/Collate.smk"
if 'catalogue' in MULTI_STEPS:
    include:
        "workflow/rules/modules/multi/Catalogue.smk"
if 'dereplicate' in MULTI_STEPS:
    include:
        "workflow/rules/modules/multi/Dereplicate.smk"
if 'DB' in MULTI_STEPS:
    include:
        "workflow/rules/modules/multi/DB.smk"
if 'visualize' in MULTI_STEPS:
    include:
        "workflow/rules/modules/multi/Vis.smk"

localrules: SamplesPrint, ALL
# master command
rule ALL:
    input:
        expand("status/{step}.done",step=MULTI_STEPS),
        "sample_table.tsv"
    output:
        touch('status/multi.done')

rule SamplesPrint:
    input:
        config['sample_table']
    output:
        "sample_table.tsv"
    run:
        samples.to_csv(path_or_buf=output[0],sep="\t",index=False,index_label=False)

rule bwa_index:
    input:
        ancient("{fasta}")
    output:
        "{fasta}.amb",
        "{fasta}.bwt",
        "{fasta}.pac",
        "{fasta}.sa",
        "{fasta}.ann"
    params:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/bwa_index.{fasta}.log"
    message: "bwa_index: Indexing {wildcards.fasta} for bwa."
    shell:
        """
        bwa index {wildcards.fasta} > {log} 2>&1
        """
