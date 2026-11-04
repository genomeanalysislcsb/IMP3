#
# kraken_base.smk
#
# This file defines generic, inheritable base rules for running and
# processing kraken2 outputs (concatenating, combining, reporting, filtering).

rule kraken2_base:
    """A generic base rule for running kraken2."""
    output:
        # This will be overridden by the inheriting rule
        temp("kraken2_base.out")
    input:
        # This will be overridden by the inheriting rule
        "reads.fq.gz"
    params:
        # These params must be defined by the inheriting rule
        database = "",
        extra_opts = "", # for --paired, --report, etc.
        # This param is optional
        tmp_out = lambda wildcards, output: str(output[0]).replace(".gz", "")
    threads:
        getThreads(128) if KRAKEN_MMAP else getThreads(KRAKTHREADS, bigmem=True, high_mem_cluster=True)
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/kraken2_base.log" # Should be overridden
    message:
        "Running base kraken2 rule." # Should be overridden
    shadow: "copy-minimal"
    shell:
        """
        kraken2 --db {params.database} {KRAKEN_MMAP_OPT} \\
                --threads {threads} \\
                --confidence {config[kraken][confidence]} \\
                --minimum-hit-groups {config[kraken][min_hit_groups]} \\
                --use-names \\
                {params.extra_opts} \\
                --output {params.tmp_out} \\
                {input} > {log} 2>&1

        # Gzip the output
        gzip {params.tmp_out} >> {log} 2>&1
        """

rule concat_kraken_base:
    """Base rule to concatenate two gzipped kraken output files."""
    input:
        pe = "pe.in.gz",
        se = "se.in.gz"
    output:
        temp("all.out.gz")
    log:
        "logs/concat_kraken_base.log"
    message:
        "Base rule: Concatenating kraken results."
    shell:
        "zcat {input.pe} {input.se} | gzip -c > {output} 2> {log}"

rule combine_kraken_base:
    """Base rule to combine GTDB and NT kraken results."""
    output:
        temp("combined.out.gz")
    input:
        gtdb = "gtdb.in.gz",
        nt = "nt.in.gz"
    params:
        script = os.path.join(SRCDIR, "combine_kraken_outputs.py"),
        unified_tax_path = config["kraken_dual_db"]["unified_tax_path"]
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/combine_kraken_base.log"
    message:
        "Base rule: Combining kraken results."
    shell:
        """
        python {params.script} \\
            --gtdb_input {input.gtdb} \\
            --nt_input {input.nt} \\
            --output {output} \\
            --unified_taxonomy {params.unified_tax_path} > {log} 2>&1
        """

rule make_report_base:
    """Base rule to generate a kraken report from a kraken output file."""
    input:
        "combined.in.gz"
    output:
        "unfiltered.report"
    params:
        database = os.path.join(DBPATH, config["kraken_dual_db"]["nt_name"]),
        report_minimizer = '--report-minimizer-data' if config['kraken']['report_minimizer'] else ''
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/make_report_base.log"
    message:
        "Base rule: Generating combined kraken report."
    shell:
        """
        kraken2 --db {params.database} {KRAKEN_MMAP_OPT} \\
                --report {output} \\
                {params.report_minimizer} \\
                {input} > {log} 2>&1
        """

rule filter_report_base:
    """Base rule to filter a kraken report."""
    input:
        "unfiltered.report"
    output:
        "filtered.report"
    params:
        min_reads = int(config['kraken']['report_min_reads']),
        min_uniq_mini = int(config['kraken']['report_min_unique_minimizers'])
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/filter_report_base.log"
    message:
        "Base rule: Filtering combined Kraken report."
    shell:
        """
        python {SRCDIR}/filter_kraken_report.py {input} \\
                                                {output} \\
                                                --min_reads {params.min_reads} \\
                                                --min_uniq_mini {params.min_uniq_mini} > {log} 2>&1
        """
