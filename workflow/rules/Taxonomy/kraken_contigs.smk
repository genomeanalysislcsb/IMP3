from os.path import getsize, join
KRAKSIZE = float(os.path.getsize(DBPATH + "/" + config['kraken']['krakendb'] + "/hash.k2d" )) / 1000000000
KRAKTHREADS = max(1, round((KRAKSIZE / (float(HIGH_MEM_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]])) * 1.15))

rule kraken_contigs:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS,
    output:
        "Analysis/taxonomy/kraken/%s.kraken.contig.output.gz" % ASS,
        "Analysis/taxonomy/kraken/%s.kraken.contig.unfiltered.report" % ASS
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


rule filter_kraken_contig_report:
    input:
        report_unfilt = "Analysis/taxonomy/kraken/%s.kraken.contig.unfiltered.report" % ASS
    output:
        report = report("Analysis/taxonomy/kraken/%s.kraken.contig.report" % ASS, category="Taxonomy")
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        script_path = os.path.join(ROOTDIR, "submodules", "KrakenTools", "combine_kreports.py"),
        min_reads = int(config['kraken']['report_min_reads']),
        min_uniq_mini = int(config['kraken']['report_min_unique_minimizers'])
    log:
        "logs/analysis/taxonomy/filter_kraken_contigs_report.log"
    message:
        "Filtering Kraken reports for contigs into a single report."
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    threads:
        getThreads(12, bigmem=False, high_mem_cluster=False)
    shell:
        """
        python {SRCDIR}/filter_kraken_report.py {input.report_unfilt} \
                                                {output.report} \
                                                --min_reads {params.min_reads} \
                                                --min_uniq_mini {params.min_uniq_mini} > {log} 2>&1
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


