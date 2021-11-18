##################################################
# Adapters
##################################################

# Adapters from Trimmomatic repository
rule trimmomatic_adapters:
    output:
        done=os.path.join("adapters", "trimmomatic.done"),
        repo=temp(directory("downloads/git/trimmomatic")),
    log:
        os.path.join("adapters", "trimmomatic.log")
    params:
        url=config["trimmomatic"]["url"],
        commit=config["trimmomatic"]["commit"],
    message:
        "Adapters: Trimmomatic adapters"
    shell:
        "(old_pwd=$(realpath ${{PWD}}) && "
        "git clone {params.url} {output.repo} && "
        "cd {output.repo} && git checkout {params.commit} && cd ${{old_pwd}} && "
        "rsync -rtP {output.repo}/adapters/ $(dirname {output.done})/ && touch {output.done}) &> $(realpath {log})"

##################################################
# rRNA
##################################################

# rRNA references from SortMeRNA repository
rule sortmerna_rrna:
    output:
        done=os.path.join("sortmerna", "sortmerna.done"),
        repo=temp(directory("downloads/git/sortmerna")),
    log:
        os.path.join("sortmerna", "sortmerna.log")
    params:
        url=config["sortmerna"]["url"],
        commit=config["sortmerna"]["commit"],
    message:
        "Filtering: SortMeRNA files"
    shell:
        "(old_pwd=$(realpath ${{PWD}}) && "
        "git clone {params.url} {output.repo} && "
        "cd {output.repo} && git checkout {params.commit} && cd ${{old_pwd}} && "
        "rsync -rtP {output.repo}/data/rRNA_databases/*.fasta $(dirname {output.done})/ && touch {output.done}) &> $(realpath {log})"

##################################################
# Genomes
##################################################

# Download given genome references
rule genome_download:
    output:
        temp("downloads/files/{genome}.{ext}"),
    log:
        "downloads/files/{genome}.{ext}.download.log"
    wildcard_constraints:
        genome="|".join(config["genomes"].keys())
    params:
        url=smk_get_genome_url
    message:
        "Download: {params.url}"
    shell:
        "wget -O {output} -o {log} --no-verbose {params.url}"

rule genome_prep:
    input:
        smk_get_genome_download_target
    output:
        temp("downloads/prep/{genome}.{ext}"),
    wildcard_constraints:
        genome="|".join(config["genomes"].keys())
    message:
        "Prepare: {input}"
    shell:
        # NOTE: if multiple extensions are present only the last one is checked
        "ifile={input} && ext=\"${{ifile##*.}}\" && "
        "if [ \"${{ext}}\" == \"gz\" ] || [ \"${{ext}}\" == \"gzip\" ]; then zcat {input} > {output}; "
        "else rsync -rtP {input} {output}; fi"

rule genome_file:
    input:
        "downloads/prep/{genome}.{ext}"
    output:
        "filtering/{genome}.{ext}"
    wildcard_constraints:
        genome="|".join(config["genomes"].keys())
    message:
        "Copy: {input}"
    shell:
        "rsync -rtP {input} {output}"