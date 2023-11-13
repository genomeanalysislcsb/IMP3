import math

localrules: allmotu

rule allmotu:
    input:
        expand("Analysis/taxonomy/mOTUs/{type}.mOTU.counts.tsv", type=TYPES)
    output:
        touch("status/taxonomy_mOTUs.done")
    threads: 1
    message:
        "allmotu: mOTUs done?"

# https://github.com/motu-tool/mOTUs/wiki
rule motus:
    input:
        r1="Preprocessing/{type}.r1.preprocessed.fq",
        r2="Preprocessing/{type}.r2.preprocessed.fq",
        se="Preprocessing/{type}.se.preprocessed.fq"
    output:
        report("Analysis/taxonomy/mOTUs/{type}.mOTU.counts.tsv",category="Taxonomy")
    threads:
        # Ensure at least 30gb of mem
        max(int(math.ceil(30 / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
    params:
        # https://github.com/motu-tool/mOTUs/wiki/Parameters-to-change-the-resulting-profiles
        # -A: print all tax. levels together
        # -c: return counts
        params = "-A -c",
        sample = config["sample"]
        # database = os.path.join(DBPATH, 'motus')
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_motus.yaml")
    log:
        "logs/analysis_motus.{type}.log"
    message:
        "motus: Running mOTUs on {wildcards.type} reads."
    shell:
        """
        if [ -s {input.se} ] && [ -s {input.r1} ] && [ -s {input.r2} ]; then
            motus profile -f {input.r1} \
                          -r {input.r2} \
                          -s {input.se} \
                          -o {output} \
                          -n {params.sample} \
                          {params.params} \
                          -t {threads} >& {log}
        elif [ -s {input.se} ]; then
            motus profile -s {input.se} \
                          -o {output} \
                          -n {params.sample} \
                          {params.params} \
                          -t {threads} >& {log}
        else
            motus profile -f {input.r1} \
                          -r {input.r2} \
                          -o {output} \
                          -n {params.sample} \
                          {params.params} \
                          -t {threads} >& {log}
        fi 
        """
