localrules: ALL_stats 

gathered = [report("Stats/all_stats.tsv",category="Stats"),
            "Stats/all_stats.Rdata"]
if "binning" in IMP_STEPS:
    gathered.append(report("Stats/all_bin_stats.tsv",category="Binning"))

rule gather_stats:
    input:
        wf_control="status/workflow.done",
        prep_control="status/preprocessing.analysed" if not "preprocessing" in IMP_STEPS else [],
        cont_control="status/assembly.analysed" if CONTIGS else []
    output:
        gathered
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        steps = IMP_STEPS,
        ass = ASS,
        types = TYPES,
        ann = config["annotation"],
        hmmer_dbs = config["hmm_DBs"] if config['annotation'] == "hmmer" else "",
        mantis_dir = "Analysis/annotation/" if config["annotation"] == "mantis" else [],
        mantis_pattern = ".mantis.+tsv.summary" if config["annotation"] == "mantis" else []
    log: "logs/stats_gather_stats.log"
    message: "gather_stats: Summarizing stats."
    threads: 1
    conda: ENVDIR + "/IMP_binning.yaml"
    script:
        SRCDIR + "/summarize_all_stats.R"
         
rule ALL_stats:
    input:
        gathered
    output:
        touch("status/stats.done")

