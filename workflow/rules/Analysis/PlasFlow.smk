##################################################
# PlasFlow is disabled

localrules: plasflow_disabled
rule plasflow_disabled:
    output:
        "Analysis/annotation/plasflow.notdone"
    message:
        "Do not run PlasFlow"
    shell:
        "touch {output}"

##################################################
# PlasFlow is enabled

# Filter input FASTA
rule plasflow_input:
    input:
        "Analysis/annotation/prokka.fna"
    output:
        temp("Analysis/annotation/prokka.plasflow.fna")
    log:
        "logs/annotation_plasflow_filter.log"
    params:
        runtime = "00:10:00",
        mem = MEMCORE,
        minlen=config["plasflow"]["minlen"]
    message:
        "Plasmid prediction w/ PlasFlow: {input}"
    shell:
        # script from PathoFact
        "{config[pathofact][scripts]}/filter.pl {params.minlen} {input} > {output} 2> {log}"

# PlasFlow output
if config["pathofact"]["enabled"]:
    # PlasFlow output from PathoFact's results
    rule plasflow_from_pathofact:
        input:
            plasflow="Analysis/annotation/pathofact/MGE/plasmid/pathofact_plasflow_aggregated.tsv",
            idmap="Analysis/annotation/pathofact/renamed/pathofact_Contig_translation.tsv"
        output:
            "Analysis/annotation/plasflow.tsv"
        params:
            runtime = "00:10:00",
            mem = MEMCORE
        threads: 1
        message:
            "Plasmid prediction w/ PlasFlow from PathoFact results: {input}"
        run:
            import re
            import pandas
            # contig ID mapping
            df_id = pandas.read_csv(input.idmap, sep="\t", header=None, names=["new", "old"], dtype=str)
            df_id.old = df_id.old.apply(lambda x: re.sub("^>", "", x))
            # PlasFlow results
            df_pf = pandas.read_csv(input.plasflow, sep="\t", header=0, index_col=0, dtype=str)
            # merge by contig ID
            df_pf = df_id.merge(
                right=df_pf,
                how="right",
                left_on="new",
                right_on="contig_name",
                validate="one_to_one"
            )
            assert all(pandas.notnull(df_pf.old))
            # select and rename columns
            df_pf = df_pf[["old", "contig_length", "label"] + [c for c in df_pf.columns if re.match("(chromosome|plasmid)\.", c)]]
            df_pf.rename(columns={"old": "contig_name"}, inplace=True)
            # save
            df_pf.to_csv(output[0], sep="\t", na_rep="NA", header=True, index=False)
else:
    # Run PlasFlow
    rule plasflow:
        input:
            "Analysis/annotation/prokka.plasflow.fna"
        output:
            "Analysis/annotation/plasflow.tsv"
        log:
            "logs/annotation_plasflow.log"
        params:
            runtime = "02:00:00",
            mem = BIGMEMCORE,
            threshold = config["plasflow"]["threshold"]
        conda:
            ENVDIR + "/PlasFlow.yaml"
        message:
            "Plasmid prediction w/ PlasFlow: {input}"
        shell:
            "PlasFlow.py --input {input} --output {output}.tmp --threshold {params.threshold} &> {log} && "
            "cut -f3,4,6- {output}.tmp > {output} && "
            "rm {output}.tmp*"
