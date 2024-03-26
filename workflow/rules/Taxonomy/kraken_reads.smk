from os.path import getsize, join
KRAKSIZE = float(os.path.getsize(DBPATH + "/" + config['kraken']['krakendb'] + "/hash.k2d" )) / 1000000000
KRAKTHREADS = max(1, round((KRAKSIZE / (float(HIGH_MEM_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]])) * 1.15))

rule kraken_reads:
    input:
        r1 = "Preprocessing/{type}.r1.preprocessed.fq",
        r2 = "Preprocessing/{type}.r2.preprocessed.fq",
    output:
        "Analysis/taxonomy/kraken/{type}.pe.kraken.output.gz",
        report("Analysis/taxonomy/kraken/{type}.pe.kraken.report", category="Taxonomy")
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = HIGH_MEM_BIGMEMCORE_MB
    params:
        tmp_out = lambda wildcards: f"Analysis/taxonomy/kraken/{wildcards.type}.pe.kraken.output",
        confidence = config['kraken']['confidence'],
        database = os.path.join(DBPATH, config['kraken']['krakendb']),
        report_minimizer = '--report-minimizer-data' if config['kraken']['report_minimizer'] else '',
        min_hit_groups = config['kraken']['min_hit_groups'],
        rrna_r1 = "Preprocessing/mt.r1.trimmed.rrna.fq.gz",
        rrna_r2 = "Preprocessing/mt.r2.trimmed.rrna.fq.gz",
        keep_rrnas = KEEP_RRNAs
    threads:
        getThreads(KRAKTHREADS, bigmem=True, high_mem_cluster=True)
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/analysis/taxonomy/kraken_reads.{type}.pe.log"
    message:
        "kraken_reads: Running kraken2 on {wildcards.type} reads against {params.database}."
    shell:
        """
        # Run kraken2 on paired-end reads
        mkdir -p Analysis/taxonomy/kraken >> {log} 2>&1

        # If `type` is `mt` add rRNAs to read files
        if [ "{wildcards.type}" == "mt" ] && [ "{params.keep_rrnas}" == "False" ]; then
            k_r1="{input.r1}.tmp"
            k_r2="{input.r2}.tmp"
            cat {input.r1} <(zcat {params.rrna_r1}) > ${{k_r1}}
            cat {input.r2} <(zcat {params.rrna_r2}) > ${{k_r2}}
        else
            k_r1="{inputr1}"
            k_r2="{input.r2}"
        fi

        kraken2 --db {params.database} \
                --confidence {params.confidence} \
                --threads {threads} \
                --paired \
                --use-names \
                --minimum-hit-groups {params.min_hit_groups} \
                {params.report_minimizer} \
                --output {params.tmp_out} \
                --report {output[1]} \
                ${{k_r1}} ${{k_r2}} >> {log} 2>&1

        # If `type` is `mt` remove temporary read files
        if [ "{wildcards.type}" == "mt" ]; then
            rm -f {input.r1}.tmp {input.r2}.tmp
        fi

        gzip {params.tmp_out} >> {log} 2>&1
        """

rule combine_kraken_reports:
    input:
        pe_report = "Analysis/taxonomy/kraken/{type}.pe.kraken.report",
        se_report = "Analysis/taxonomy/kraken/{type}.se.kraken.report"
    output:
        report("Analysis/taxonomy/kraken/{type}.kraken.report",category="Taxonomy")
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        script_path = os.path.join(ROOTDIR, "submodules", "KrakenTools", "combine_kreports.py"),
        combined_report_options = "--only-combined --no-headers",
        has_minimizer = '--has-minimizer' if config['kraken']['report_minimizer'] else ''
    log:
        "logs/analysis/taxonomy/combine_kraken_reports.{type}.log"
    message:
        "Combining Kraken reports for {wildcards.type} into a single report."
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    threads:
        getThreads(12, bigmem=False, high_mem_cluster=False)
    shell:
        """
        python {params.script_path} \
            -r {input.pe_report} {input.se_report} \
            -o {output} \
            {params.combined_report_options} \
            {params.has_minimizer} > {log} 2>&1
        """

# Run SeqKit to get read stats
rule seqkit_stats:
    input:
        "Preprocessing/{type}.r1.preprocessed.fq",
        "Preprocessing/{type}.r2.preprocessed.fq",
        "Preprocessing/{type}.se.preprocessed.fq"
    output:
        "Stats/{type}/{type}.stats"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log:
        "logs/analysis/taxonomy/seqkit_stats.{type}.log"
    message:
        "seqkit_stats: Running seqkit stats on {wildcards.type} reads."
    conda:
        os.path.join(ENVDIR, "IMP_assembly.yaml")
    threads:
        getThreads(12, bigmem=False, high_mem_cluster=False)
    shell:
        """
        seqkit stats --threads {threads} \
                     --all \
                     --tabular \
                     --basename \
                     {input} > {output} 2>> {log}
        """

# Helper function to read the seqkit stats file, check which kmer_distrib file is the closest to the min read length
def get_kmer_distrib(ome):
    import os
    import glob

    # Get the min length from the seqkit stats file
    min_len = None
    with open(f"Stats/{ome}/{ome}.stats", "r") as f:
        header = next(f)
        for line in f:
            line = line.strip().split("\t")
            # Min length is in column 6
            if min_len is None or int(line[5]) < min_len:
                min_len = int(line[5])
    
    # Get the kmer_distrib file that is closest to the min length
    kmer_distrib_files = glob.glob(os.path.join(DBPATH, config['kraken']['krakendb'], "*.kmer_distrib"))

    if not kmer_distrib_files:
        print("Warning: No kmer_distrib files found in the database directory.")
        return None
    
    diff = None
    closest_kmer_distrib_file_size = None
    for kmer_distrib_file in kmer_distrib_files:
        kmer_distrib_file_name = os.path.basename(kmer_distrib_file)
        kmer_distrib_file_size = int(kmer_distrib_file_name.split(".")[0].replace("database", "").replace("mers", ""))
        diff_current = abs(kmer_distrib_file_size - min_len)
        if closest_kmer_distrib_file_size is None:
            closest_kmer_distrib_file_size = kmer_distrib_file_size
            diff = diff_current
        elif diff_current < diff:
            closest_kmer_distrib_file_size = kmer_distrib_file_size
            diff = diff_current
    
    if diff > 25:
        print(f"Warning: The closest kmer_distrib file to the min read length is more than 25bp away with a difference of {diff}bp.")
    
    return closest_kmer_distrib_file_size

# TODO: Improve stability of this, eg with checkpoint rule to properly handle missing kmer_distrib files
rule bracken_reads:
    input:
        seqkit_stats = "Stats/{type}/{type}.stats",
        kraken_report ="Analysis/taxonomy/kraken/{type}.kraken.report"
    output:
        bracken_dir = directory("Analysis/taxonomy/bracken/{type}"),
        check = "status/{type}.bracken.done"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        database = os.path.join(DBPATH, config['kraken']['krakendb']),
        kmer_distrib_file = lambda wildcards: get_kmer_distrib(wildcards.type),
        level = config['bracken']['levels'],
        min_reads = config['bracken']['min_reads'],
        ome = lambda wildcards: wildcards.type,
        kraken_dir = "Analysis/taxonomy/kraken"
    log:
        "logs/analysis/taxonomy/bracken_reads.{type}.log"
    message:
        "bracken_reads: Running bracken on {wildcards.type} reads."
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    threads:
        getThreads(32, bigmem=False, high_mem_cluster=False)
    shell:
        """
        # Ensure Bracken dir exists
        mkdir -p {output.bracken_dir} >> {log}

        # Just touch output file if params.kmer_distrib is None (i.e. no kmer_distrib file found)
        if [ -z "{params.kmer_distrib_file}" ]; then
            echo "Warning: No kmer_distrib file found for {wildcards.type}. Skipping bracken." >> {log}
            touch "status/{params.ome}.bracken.not_run"
            touch {output.check}
        else
            # Iterate over all levels
            for level in {params.level}; do
                # Run bracken on kraken report, if it fails touch output and continue
                bracken -d {params.database} \
                        -i {input.kraken_report} \
                        -o {output.bracken_dir}/{params.ome}.${{level}}.bracken \
                        -r {params.kmer_distrib_file} \
                        -l ${{level}} \
                        -t {threads} >> {log} 2>&1 || touch {output.check} && continue
                
                # Check if bracken output exists and isn't empty, if it isnt move the report from kraken dir to bracken dir
                if [ ! -s "{output.bracken_dir}/{params.ome}.${{level}}.bracken" ]; then
                    echo "Error: Bracken output for {wildcards.type} at level ${{level}} is empty." >> {log}
                    exit 1
                else
                    mv {params.kraken_dir}/{params.ome}.kraken_bracken_*.report {output.bracken_dir}/
                fi
            done

            touch {output.check}
        fi
        echo "Finished Bracken for {wildcards.type}." >> {log}
        """
