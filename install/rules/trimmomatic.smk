##################################################
# Trimmomatic adapters
##################################################

rule trimmomatic_adapters:
    input:
        "downloads/git/trimmomatic"
    output:
        os.path.join("adapters", "trimmomatic.done")
    log:
        os.path.join("adapters", "trimmomatic.log")
    message:
        "Adapters: Trimmomatic adapters"
    shell:
        "rsync -rtP --log-file={log} {input}/adapters/ $(dirname {output})/ && touch {output}"
