##################################################
# SortMeRNA rRNA references
##################################################

rule sortmerna_rrna:
    input:
        "downloads/git/sortmerna"
    output:
        os.path.join("sortmerna", "sortmerna.done")
    log:
        os.path.join("sortmerna", "sortmerna.log")
    message:
        "Filtering: SortMeRNA files"
    shell:
        "rsync -rtP --log-file={log} {input}/data/rRNA_databases/*.fasta $(dirname {output})/ && touch {output}"
