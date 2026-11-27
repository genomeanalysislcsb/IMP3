##################################################
# Deacon human pangenome index (host read removal)
##################################################
rule deacon_index_fetch:
    output:
        "filtering/deacon/{index}.k31w15.idx"
    wildcard_constraints:
        index = config["deacon"]["index"]
    log:
        "filtering/deacon/{index}.log"
    conda:
        os.path.join(ENV_DIR, "deacon.yaml")
    message:
        "Filtering: Deacon index {wildcards.index}"
    shell:
        "cd $(dirname {output}) && deacon index fetch {wildcards.index} > $(basename {log}) 2>&1 && "
        "test -s $(basename {output})"

##################################################
# Deacon indexes of the other filter references (PhiX, ...), filtering: method: deacon
##################################################
rule deacon_index_build:
    input:
        "filtering/{genome}.fa"
    output:
        "filtering/deacon/{genome}.k31w15.idx"
    wildcard_constraints:
        genome = "|".join(config.get("genomes", {}).keys()) or "no_genome"
    log:
        "filtering/deacon/{genome}.log"
    conda:
        os.path.join(ENV_DIR, "deacon.yaml")
    message:
        "Filtering: Deacon index of {wildcards.genome}"
    shell:
        "deacon index build {input} > {output}.tmp 2> {log} && mv {output}.tmp {output}"
