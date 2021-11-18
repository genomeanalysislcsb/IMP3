##################################################
# Kraken2 databases
##################################################
# NOTE: split into separate rules
#   - to include other decompression rules if needed
#   - archives may contain a folder structure (which should be ignored)

rule kraken2db_download:
    output:
        temp("downloads/files/{kdb}.{ext}")
    log:
        "downloads/files/{kdb}.{ext}.log"
    wildcard_constraints:
        kdb="|".join(config["kraken2"].keys())
    params:
        url=smk_get_kraken2db_url
    threads: 1
    message:
        "Kraken2 database: download {output} from {params.url}"
    shell:
        "wget -O {output} -o {log} --no-verbose {params.url}"

rule kraken2db_decompression:
    input:
        smk_get_kraken2db_download_target
    output:
        temp(directory("downloads/decompressed/{kdb}"))
    log:
        "downloads/decompressed/{kdb}.log"
    wildcard_constraints:
        kdb="|".join(config["kraken2"].keys())
    message:
        "Kraken2 database: decompress {input} into {output}"
    shell:
        # NOTE: if multiple extensions are present only the last one is checked
        "ifile={input} && ext=\"${{ifile##*.}}\" && "
        "mkdir -p {output} && "
        "if [ \"${{ext}}\" == \"tgz\" ]; then tar -xzvf {input} -C {output} &> {log}; fi"

rule kraken2db_files:
    input: # folder, decompressed output
        "downloads/decompressed/{kdb}"
    output:
        files=expand("{{kdb}}/{fname}", fname=get_kraken2db_file_names()),
        flist=temp("{kdb}/files.txt"),
    log:
        "{kdb}/download.log"
    wildcard_constraints:
        kdb="|".join(config["kraken2"].keys())
    message:
        "Kraken DB: {wildcards.kdb}: {output}"
    shell:
        "find {input}/ -type f ! -name '.snakemake_timestamp' | sed 's@{input}/@@' | sort > {output.flist} && "
        "rsync -rtP --no-R --log-file={log} --files-from={output.flist} {input} $(dirname {output.flist})/ 2>> {log}"

##################################################
# GTDB-tk
##################################################

rule download_gtdbtk:
    output:
        temp("downloads/files/gtdbtk.{ext}".format(ext=config["gtdbtk"]["compression"]))
    log:
        "gtdbtk/download.log"
    params:
        url=config["gtdbtk"]["url"]
    message:
        "Download GTDB-tk data: {output}"
    shell:
        "wget -O {output} -o {log} --no-verbose {params.url} && "

rule decompress_gtdbtk:
    input:
        temp("downloads/files/gtdbtk.{ext}".format(ext=config["gtdbtk"]["compression"]))
    output:
        "gtdbtk/download.done"
    log:
        "gtdbtk/decompression.log"
    params:
        url=config["gtdbtk"]["url"]
    message:
        "Decompress GTDB-tk data: {output}"
    shell:
        "tar -xzvf {input} -C $(dirname {output}) $> {log} && touch {output}"
