import math

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
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        steps = IMP_STEPS,
        ass = ASS,
        types = TYPES,
        ann = config["annotation"],
        hmmer_dbs = config["hmm_DBs"] if config['annotation'] == "hmmer" else "",
        mantis_dir = "Analysis/annotation/" if config["annotation"] == "mantis" else [],
        mantis_pattern = ".mantis.+tsv.summary" if config["annotation"] == "mantis" else []
    log:
        "logs/misc/stats_gather_stats.log"
    message:
        "gather_stats: Summarizing stats."
    threads:
        # Ensure at least 5gb of mem
        max(int(math.ceil(5 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    script:
        SRCDIR + "/summarize_all_stats.R"
         
rule ALL_stats:
    input:
        gathered
    output:
        touch("status/stats.done")

