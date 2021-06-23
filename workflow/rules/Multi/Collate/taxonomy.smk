checkpoint collate_kraken:
    input:
        lambda wildcards: expand("{raw_directory}/Analysis/taxonomy/kraken/{{type}}.kraken.report",raw_directory=samples.loc[[a and b for a,b in zip(samples.taxonomy == 1,samples[wildcards.type] == 1)],"path"])
    output:
        "Collection/multi_{type}.kraken.Rdata"
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE 
    params:
        samples = lambda wildcards: samples.loc[[a and b for a,b in zip(samples.taxonomy == 1,samples[wildcards.type] == 1)],"sample"],
        out_prefix = lambda wildcards: "Collection/multi_" + wildcards.type + ".kraken."
    log: "logs/collate_kraken_{type}.kraken.log"
    message: 'collate_kraken: collecting kraken {wildcards.type} data.'
    threads: 1
    conda: ENVDIR + "/IMP_multi.yaml"
    script:
        SRCDIR + "/multi_collect_kraken.R"

rule collate_mOTUs:
    input:
        lambda wildcards: expand("{raw_directory}/Analysis/taxonomy/mOTUs/{{type}}.mOTU.counts.tsv",raw_directory=samples.loc[[a and b for a,b in zip(samples.taxonomy == 1,samples[wildcards.type] == 1)],"path"])
    output:
        "Collection/multi_{type}.mOTUs.Rdata",
        expand("Collection/multi_{{type}}.mOTUs.{rank}.tsv",rank=["k","p","c","o","f","g","s"])
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        samples = lambda wildcards: samples.loc[[a and b for a,b in zip(samples.taxonomy == 1,samples[wildcards.type] == 1)],"sample"]
    log: "logs/collate_mOTUs_{type}.mOTUs.log"
    message: 'collate_mOTUs: collecting mOTUs {wildcards.type} data.'
    threads: 1
    conda: ENVDIR + "/IMP_multi.yaml"
    script:
        SRCDIR + "/multi_collect_mOTUs.R"
        

checkpoint collate_EukDetect:
    input:
        lambda wildcards: expand("{raw_directory}/Analysis/taxonomy/EukDetect/{{type}}_filtered_hits_table.txt",raw_directory=samples.loc[[a and b for a,b in zip(samples.EukDetect == 1,samples[wildcards.type] == 1)],"path"])
    output:
        "Collection/multi_{type}.EukDetect.Rdata",
        expand("Collection/multi_{{type}}.EukDetect.{kind}.tsv",kind=["Read_counts","Total_marker_coverage","Percent_identity"])
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        samples = lambda wildcards: samples.loc[[a and b for a,b in zip(samples.EukDetect == 1,samples[wildcards.type] == 1)],"sample"],
        out_prefix = lambda wildcards: "Collection/multi_" + wildcards.type + ".EukDetect."
    log: "logs/collate_EukDetect_{type}.EukDetect.log"
    message: 'collate_EukDetect: collecting EukDetect {wildcards.type} data.'
    threads: 1
    conda: ENVDIR + "/IMP_multi.yaml"
    script:
        SRCDIR + "/multi_collect_eukdetect.R"

localrules:  ctrl_tax_collation

classifiers=["kraken","mOTUs"]
if 1 in samples.EukDetect.unique():
    classifiers.append("EukDetect")

rule ctrl_tax_collation:
    input:
        expand("Collection/multi_{type}.{classifier}.Rdata",type=TYPES,classifier=classifiers)
    output:
        touch("status/taxonomy.collated")

