import math

localrules: allmotu_se

rule allmotu_se:
    input:
        expand("Analysis/taxonomy/mOTUs/{type}.mOTU.{out_type}.tsv", type=TYPES, out_type=["counts", "rel_abun"])
    output:
        touch("status/taxonomy_mOTUs.done")
    threads:
        1
    message:
        "allmotu: mOTUs done?"

rule get_motus_db:
    output:
        os.path.join(ENVDIR, 'motus.ready')
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        4
    conda:
        os.path.join(ENVDIR, "IMP_motus.yaml")
    log:
        os.path.join(ENVDIR, "logs/motus_downloadDB.log")
    message:
        "motus: Downloading mOTUs database."
    shadow: "copy-minimal"
    shell:
        """
        motus downloadDB > {log} 2>&1 \
         && touch {output} >> {log} 2>&1
        """

def set_motu_runtime(wildcards, attempt):
    if attempt == 1:
        return VERY_SHORT_RUNTIME
    elif attempt == 2:
        return SHORTISH_RUNTIME
    else:
        return MEDIUM_RUNTIME

# https://github.com/motu-tool/mOTUs/wiki
rule motus:
    input:
        db_ready = os.path.join(ENVDIR, 'motus.ready'),
        se="Preprocessing/{type}.se.preprocessed.fq.gz"
    output:
        report("Analysis/taxonomy/mOTUs/{type}.mOTU.counts.tsv", category="Taxonomy")
    threads:
        # Ensure at least 150gb of mem
        max(int(math.ceil(150 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
    params:
        # https://github.com/motu-tool/mOTUs/wiki/Parameters-to-change-the-resulting-profiles
        # -A: print all tax. levels together
        # -c: return counts
        extra_params = "-A -c",
        min_markers = config["motus"]["min_markers"],
        sample = config["sample"],
        min_size = lambda wildcards: config["motus"][f"min_read_len_{wildcards.type}"]
        # database = os.path.join(DBPATH, 'motus')
    resources:
        runtime = set_motu_runtime,
        mem_mb = MAIN_MEMCORE_MB
    retries:
        3
    conda:
        os.path.join(ENVDIR, "IMP_motus.yaml")
    log:
        "logs/analysis/taxonomy/motus.{type}.counts.log"
    message:
        "motus: Running mOTUs (counts) on {wildcards.type} reads."
    shadow: "copy-minimal"
    shell:
        """
        motus profile -s {input.se} \
                      -o {output} \
                      -n {params.sample} \
                      {params.extra_params} \
                      -g {params.min_markers} \
                      -l {params.min_size} \
                      -t {threads} >& {log}
        """

use rule motus as motus_rel_ab with:
    output:
        report("Analysis/taxonomy/mOTUs/{type}.mOTU.rel_abun.tsv",category="Taxonomy")
    params:
        extra_params = "-p -A",
        min_markers = config["motus"]["min_markers"],
        sample = config["sample"],
        min_size = lambda wildcards: config["motus"][f"min_read_len_{wildcards.type}"]
    log:
        "logs/analysis/taxonomy/motus.{type}.rel_abun.log"
    message:
        "motus: Running mOTUs (relative abundance) on {wildcards.type} reads."