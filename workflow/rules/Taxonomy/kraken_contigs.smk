from os.path import getsize, join
KRAKSIZE = float(os.path.getsize(DBPATH + "/" + config['kraken']['krakendb'] + "/hash.k2d" )) / 1000000000
KRAKTHREADS = max(1, round((KRAKSIZE / (float(HIGH_MEM_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]])) * 1.15))

rule kraken_contigs:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS,
    output:
        "Analysis/taxonomy/kraken/%s.kraken.contig.output.gz" % ASS,
        "Analysis/taxonomy/kraken/%s.kraken.contig.report" % ASS
    resources:
        # runtime = SHORT_RUNTIME,
        # TODO: Change back after prospectomics done
        runtime = 360,
        mem_mb = HIGH_MEM_BIGMEMCORE_MB
    params:
        tmp_out = lambda wildcards: f"Analysis/taxonomy/kraken/{ASS}.kraken.contig.output",
        confidence = config['kraken']['confidence'],
        database = os.path.join(DBPATH, config['kraken']['krakendb']),
        report_minimizer = '--report-minimizer-data' if config['kraken']['report_minimizer'] else '',
        min_hit_groups = config['kraken']['min_hit_groups']
    threads:
        getThreads(KRAKTHREADS, bigmem=True, high_mem_cluster=True)
    conda:
        ENVDIR + "/IMP_taxonomy.yaml"
    log:
        "logs/analysis/taxonomy/kraken_contigs.log"
    message:
        "kraken_contigs: Running kraken2 on contigs against {params.database}."
    shell:
        """
        kraken2 --db {params.database} \
                --confidence {params.confidence} \
                --threads {threads} \
                --use-names \
                --minimum-hit-groups {params.min_hit_groups} \
                {params.report_minimizer} \
                --output {params.tmp_out} \
                --report {output[1]} \
                {input} > {log} 2>&1
        gzip {params.tmp_out} >> {log} 2>&1
        """


# TODO: Removed for now since we filter much stricter than previously. discuss if still useful
# rule kraken_parse:
#     input:
#         "Analysis/taxonomy/kraken/%s.kraken.contig.output" % ASS,
#         "Analysis/taxonomy/kraken/%s.kraken.contig.report" % ASS
#     output:
#         "Analysis/taxonomy/kraken/%s.kraken.parsed.tsv" % ASS
#     threads:
#         1
#     resources:
#         runtime = VERY_SHORT_RUNTIME,
#         mem_mb = MAIN_MEMCORE_MB
#     conda:
#         os.path.join(ENVDIR, "IMP_taxonomy.yaml")
#     log:
#         "logs/analysis/taxonomy/kraken_parse.log"
#     message:
#         "kraken_parse: Parsing kraken2 output for contigs."
#     shell:
#         """
#         {SRCDIR}/krakenContigEntropy_dynamic.py -c {input[0]} \
#                                                 -e 0.1 \
#                                                 -o {output} \
#                                                 -r {input[1]} > {log} 2>&1
#         """


