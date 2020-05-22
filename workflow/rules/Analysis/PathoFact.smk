##################################################
# PathoFact is disabled

localrules: pathofact_disabled
rule pathofact_disabled:
    output:
        "Analysis/annotation/pathofact.notdone"
    message:
        "Do not run PathoFact"
    shell:
        "touch {output}"

##################################################
# PathoFact is enabled

# Include rules from PathoFact
if config["pathofact"]["enabled"]:
    include:
        MODULDIR + "/PathoFact/workflows/Combine_PathoFact_workflow.smk"

# Input files
# FNA: link to file from Prokka
localrules: pathofact_input_fna
rule pathofact_input_fna:
    input:
        "Analysis/annotation/prokka.fna"
    output:
        temp("Analysis/annotation/pathofact.fna")
    message:
        "PathoFact input: {output} => {input}"
    shell:
        "cd $(dirname {output}) && ln -sf $(basename {input}) $(basename {output})"

# FAA: process file from Prokka (rm FASTA record description)
rule pathofact_input_faa:
    input:
        "Analysis/annotation/prokka.faa"
    output:
        temp("Analysis/annotation/pathofact.faa")
    params:
        runtime = "00:10:00",
        mem = MEMCORE
    threads: 1
    message:
        "PathoFact input: {output} from {input}"
    shell:
        "sed '/^>/ s/ .*//' {input} > {output}"

# Contig/gene ID mapping: Create from proc. GFF
rule pathofact_input_contig:
    input:
        "Analysis/annotation/annotation.filt.gff"
    output:
        temp("Analysis/annotation/pathofact.contig")
    params:
        runtime = "00:10:00",
        mem = MEMCORE
    threads: 1
    message:
        "PathoFact input: {output} from {input}"
    shell:
        "cut -f1,9 {input} | cut -d';' -f1 | sed 's/ID=//' > {output}"

# Output files
rule pathofact_summary:
    input:
        "Analysis/annotation/pathofact/PathoFact_pathofact_predictions.csv"
    output:
        "Analysis/annotation/pathofact.tsv"
    params:
        runtime = "00:10:00",
        mem = MEMCORE
    threads: 1
    message:
        "Pathofact summary from {input}"
    shell:
        "sed 's/,/\\t/g' {input} > {output}"

# Compress PathoFact's output
rule pathofact_compress:
    input:
        ["Analysis/annotation/pathofact.tsv", "Analysis/annotation/plasflow.tsv"] if config["plasflow"]["enabled"] else "Analysis/annotation/pathofact.tsv" 
    output:
        "Analysis/annotation/pathofact.tar.gz"
    log:
        "logs/annotation_tar_pathofact.log"
    params:
        runtime = "01:00:00",
        mem = MEMCORE,
        tar_dir = "Analysis/annotation/pathofact"
    threads: 1
    message:
        "Compress PathoFact's output directory: {params.tar_dir}"
    shell:
        "tar cvzf {output} {params.tar_dir} &> {log} && rm -r {params.tar_dir} &>> {log}"
