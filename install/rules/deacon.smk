##################################################
# Deacon human pangenome index (host read removal)
##################################################
rule deacon_index_fetch:
    output:
        "filtering/deacon/{index}.k31w15.idx"
    log:
        "filtering/deacon/{index}.log"
    conda:
        os.path.join(ENV_DIR, "deacon.yaml")
    message:
        "Filtering: Deacon index {wildcards.index}"
    shell:
        "cd $(dirname {output}) && deacon index fetch {wildcards.index} > $(basename {log}) 2>&1 && "
        "test -s $(basename {output})"
