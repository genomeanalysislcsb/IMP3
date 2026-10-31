import math

localrules: allmotu

rule allmotu:
    input:
        expand("Analysis/taxonomy/mOTUs/{type}.mOTU.{out_type}.tsv", type=TYPES, out_type=["counts", "rel_abun"])
    output:
        touch("status/taxonomy_mOTUs.done")
    threads:
        1
    message:
        "allmotu: mOTUs done?"

rule install_motus4:
    output:
        motus4_path = directory(os.path.join(SUBMOD_DIR, "motus-tool")), 
        motus_tool = os.path.join(SUBMOD_DIR, 'motus-tool/motus.py'),
        motus4_db = directory(os.path.join(SUBMOD_DIR, "motus-tool/db_mOTU")),
    resources:
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        4
    params:
        m4_src = "https://sunagawalab.ethz.ch/share/MOTUS/motus-tool/v4.0.0a/mOTUsv4.0.0a.tar.gz",
        motus4 = os.path.join(SUBMOD_DIR, 'motus-tool/motus.py'),
    conda:
        os.path.join(ENVDIR, "IMP_motus.yaml")
    log:
        os.path.join(ENVDIR, "logs/motus4_install.log")
    message:
        "motus: Installing mOTUs4."
    shadow: "copy-minimal"
    shell:
        """
        echo "Installing mOTUs4 from {params.m4_src} ..." > {log}

        cd {SUBMOD_DIR} >>& {log}
        wget {params.m4_src} >>& {log}
        tar -xzvf mOTUsv4.0.0a.tar.gz >>& {log}

        echo "Downloading mOTUs4 db ..." > {log}
        python {params.motus4} downloadDB >> {log} 2>&1

        echo "Done." >> {log}
        """

# rule get_motus_db:
#     input:
#         check = rules.install_motus4.output
#     output:
#         motus4_db = directory(os.path.join(SUBMOD_DIR, "motus-tool/db_mOTU")),
#     resources:
#         runtime = SHORT_RUNTIME,
#         mem_mb = MAIN_MEMCORE_MB
#     threads:
#         4
#     conda:
#         os.path.join(ENVDIR, "IMP_motus.yaml")
#     log:
#         os.path.join(ENVDIR, "logs/motus_downloadDB.log")
#     message:
#         "motus: Downloading mOTUs database."
#     shadow: "copy-minimal"
#     shell:
#         """
#         motus downloadDB > {log} 2>&1 \
#          && touch {output} >> {log} 2>&1
#         """

# https://github.com/motu-tool/mOTUs/wiki
rule motus:
    input:
        check_install = rules.install_motus4.output,
        r1="Preprocessing/{type}.r1.preprocessed.fq.gz",
        r2="Preprocessing/{type}.r2.preprocessed.fq.gz",
        se="Preprocessing/{type}.se.preprocessed.fq.gz"
    output:
        counts = report("Analysis/taxonomy/mOTUs/{type}.mOTU.counts.tsv", category="Taxonomy"),
        rel_ab = report("Analysis/taxonomy/mOTUs/{type}.mOTU.rel_abun.tsv",category="Taxonomy")
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
        motus4 = os.path.join(SUBMOD_DIR, 'motus-tool/motus.py'),
        base_out = lambda wildcards: "Analysis/taxonomy/mOTUs/{wildcards.type}.mOTU"
        ome = lambda wildcards: "{wildcards.type}"
    resources:
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_motus.yaml")
    log:
        "logs/analysis/taxonomy/motus.{type}.counts.log"
    message:
        "motus: Running mOTUs (counts) on {wildcards.type} reads."
    shadow: "copy-minimal"
    shell:
        """
        echo "Running mOTUs4 for ome {params.ome} ..." > {log}

        if [ -s {input.se} ] && [ -s {input.r1} ] && [ -s {input.r2} ]; then
            python {params.motus4} profile \
                -f {input.r1} \
                -r {input.r2} \
                -s {input.se} \
                -o {params.base_out} \
                -n {params.sample} \
                {params.extra_params} \
                -g {params.min_markers} \
                -l {params.min_size} \
                -t {threads} >>& {log}
        elif [ -s {input.se} ]; then
            python {params.motus4} profile \
                -s {input.se} \
                -o {params.base_out} \
                -n {params.sample} \
                {params.extra_params} \
                -g {params.min_markers} \
                -l {params.min_size} \
                -t {threads} >>& {log}
        else
            python {params.motus4} profile \
                -f {input.r1} \
                -r {input.r2} \
                -o {params.base_out} \
                -n {params.sample} \
                {params.extra_params} \
                -g {params.min_markers} \
                -l {params.min_size} \
                -t {threads} >>& {log}
        fi

        # Move outputs to final files
        echo "Moving outputs to final output files ..." >>& {log}
        mv {params.base_out} {output.counts} >>& {log}
        mv {params.base_out}.relab {output.rel_ab} >>& {log}

        echo "Done." >>& {log}
        """

# use rule motus as motus_rel_ab with:
#     output:
#         report("Analysis/taxonomy/mOTUs/{type}.mOTU.rel_abun.tsv",category="Taxonomy")
#     params:
#         extra_params = "-p -A",
#         min_markers = config["motus"]["min_markers"],
#         sample = config["sample"],
#         min_size = lambda wildcards: config["motus"][f"min_read_len_{wildcards.type}"]
#     log:
#         "logs/analysis/taxonomy/motus.{type}.rel_abun.log"
#     message:
#         "motus: Running mOTUs (relative abundance) on {wildcards.type} reads."