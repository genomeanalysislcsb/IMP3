include:
    "../../Multi/Collate/features.smk"
include:
    "../../Multi/Collate/taxonomy.smk"
include:
    "../../Multi/Collate/stats.smk"

# master command
rule all_collate:
    input:
        expand("status/{task}.collated",task=config["collection"]["results"].split())
    output:
        touch('status/collate.done')


