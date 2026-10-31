#
# kraken_reads_dual.smk
#
# This file contains the complete Snakemake rules for the experimental
# dual-database Kraken2 workflow for sequencing reads, using rule inheritance.
#

# --- Define Kraken2 runs for each read type and database ---

use rule kraken2_base as kraken_pe_gtdb with:
    input: r1="Preprocessing/{type}.r1.preprocessed.fq.gz", r2="Preprocessing/{type}.r2.preprocessed.fq.gz"
    output: temp("Analysis/taxonomy/kraken/{type}.pe.gtdb.kraken.output.gz")
    params: database=os.path.join(DBPATH, config["kraken_dual_db"]["gtdb_name"]), extra_opts="--paired"
    log: "logs/analysis/taxonomy/kraken_reads_dual.{type}.pe.gtdb.log"
    message: "Dual-DB Kraken (GTDB): Running kraken2 on PE {wildcards.type} reads."

use rule kraken2_base as kraken_pe_nt with:
    input: r1="Preprocessing/{type}.r1.preprocessed.fq.gz", r2="Preprocessing/{type}.r2.preprocessed.fq.gz"
    output: temp("Analysis/taxonomy/kraken/{type}.pe.nt.kraken.output.gz")
    params: database=os.path.join(DBPATH, config["kraken_dual_db"]["nt_name"]), extra_opts="--paired"
    log: "logs/analysis/taxonomy/kraken_reads_dual.{type}.pe.nt.log"
    message: "Dual-DB Kraken (NT): Running kraken2 on PE {wildcards.type} reads."

use rule kraken2_base as kraken_se_gtdb with:
    input: "Preprocessing/{type}.se.preprocessed.fq.gz"
    output: temp("Analysis/taxonomy/kraken/{type}.se.gtdb.kraken.output.gz")
    params: database=os.path.join(DBPATH, config["kraken_dual_db"]["gtdb_name"]), extra_opts=""
    log: "logs/analysis/taxonomy/kraken_reads_dual.{type}.se.gtdb.log"
    message: "Dual-DB Kraken (GTDB): Running kraken2 on SE {wildcards.type} reads."

use rule kraken2_base as kraken_se_nt with:
    input: "Preprocessing/{type}.se.preprocessed.fq.gz"
    output: temp("Analysis/taxonomy/kraken/{type}.se.nt.kraken.output.gz")
    params: database=os.path.join(DBPATH, config["kraken_dual_db"]["nt_name"]), extra_opts=""
    log: "logs/analysis/taxonomy/kraken_reads_dual.{type}.se.nt.log"
    message: "Dual-DB Kraken (NT): Running kraken2 on SE {wildcards.type} reads."

# --- Concatenate PE and SE results for each database ---

use rule concat_kraken_base as concat_kraken_reads_gtdb with:
    input: pe="Analysis/taxonomy/kraken/{type}.pe.gtdb.kraken.output.gz", se="Analysis/taxonomy/kraken/{type}.se.gtdb.kraken.output.gz"
    output: temp("Analysis/taxonomy/kraken/{type}.all_reads.gtdb.kraken.output.gz")
    log: "logs/analysis/taxonomy/concat_kraken_reads.{type}.gtdb.log"
    message: "Concatenating GTDB kraken results for {wildcards.type} reads."

use rule concat_kraken_base as concat_kraken_reads_nt with:
    input: pe="Analysis/taxonomy/kraken/{type}.pe.nt.kraken.output.gz", se="Analysis/taxonomy/kraken/{type}.se.nt.kraken.output.gz"
    output: temp("Analysis/taxonomy/kraken/{type}.all_reads.nt.kraken.output.gz")
    log: "logs/analysis/taxonomy/concat_kraken_reads.{type}.nt.log"
    message: "Concatenating NT kraken results for {wildcards.type} reads."

# --- Combine results and generate final report ---

use rule combine_kraken_base as combine_kraken_reads with:
    input: gtdb="Analysis/taxonomy/kraken/{type}.all_reads.gtdb.kraken.output.gz", nt="Analysis/taxonomy/kraken/{type}.all_reads.nt.kraken.output.gz"
    output: temp("Analysis/taxonomy/kraken/{type}.all_reads.combined.kraken.output.gz")
    log: "logs/analysis/taxonomy/combine_kraken_reads.{type}.log"
    message: "Combining GTDB and NT kraken results for {wildcards.type} reads."

use rule make_report_base as make_report_reads with:
    input: "Analysis/taxonomy/kraken/{type}.all_reads.combined.kraken.output.gz"
    output: "Analysis/taxonomy/kraken/{type}.kraken.unfiltered.report"
    log: "logs/analysis/taxonomy/make_report_reads.{type}.log"
    message: "Generating final combined kraken report for {wildcards.type} reads."

use rule filter_report_base as filter_kraken_read_report with:
    input: "Analysis/taxonomy/kraken/{type}.kraken.unfiltered.report"
    output: report("Analysis/taxonomy/kraken/{type}.kraken.report", category="Taxonomy")
    log: "logs/analysis/taxonomy/filter_kraken_reports.{type}.log"
    message: "Filtering combined Kraken report for {wildcards.type}."

# --- Run Bracken ---

# Run SeqKit to get read stats
rule seqkit_stats:
    input:
        r1 = "Preprocessing/{type}.r1.preprocessed.fq.gz",
        r2 = "Preprocessing/{type}.r2.preprocessed.fq.gz",
        se = "Preprocessing/{type}.se.preprocessed.fq.gz"
    output:
        "Stats/{type}/{type}.stats"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log:
        "logs/analysis/taxonomy/seqkit_stats.{type}.log"
    message:
        "seqkit_stats: Running seqkit stats on {wildcards.type} reads."
    conda:
        os.path.join(ENVDIR, "IMP_assembly.yaml")
    threads:
        getThreads(12, bigmem=False, high_mem_cluster=False)
    shadow: "minimal"
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
# TODO: bracken -r input is an int, the var is implying file, check and fix
bracken_levels_cfg = config['bracken']['levels'].split()

rule bracken_reads:
    input:
        seqkit_stats = "Stats/{type}/{type}.stats",
        kraken_report ="Analysis/taxonomy/kraken/{type}.kraken.report"
    output:
        bracken_dir = directory("Analysis/taxonomy/bracken/{type}"),
        bracken = expand("Analysis/taxonomy/bracken/{{type}}/{level}.bracken", level=bracken_levels_cfg),
        check = "status/{type}.bracken.done"
    resources:
        runtime = SHORTISH_RUNTIME,
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
    shadow: "copy-minimal"
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
                        -o {output.bracken_dir}/${{level}}.bracken \
                        -r {params.kmer_distrib_file} \
                        -l ${{level}} \
                        -t {params.min_reads} >> {log} 2>&1 || touch {output.check} && continue
                
                # Check if bracken output exists and isn't empty, if it isnt move the report from kraken dir to bracken dir
                if [ ! -s "{output.bracken_dir}/${{level}}.bracken" ]; then
                    echo "Error: Bracken output for {wildcards.type} at level ${{level}} is empty." >> {log}
                    exit 1
                else
                    mv {params.kraken_dir}/kraken_bracken_*.report {output.bracken_dir}/
                fi
            done

            touch {output.check}
        fi
        echo "Finished Bracken for {wildcards.type}." >> {log}
        """