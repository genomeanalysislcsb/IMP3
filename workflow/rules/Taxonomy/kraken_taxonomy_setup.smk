# kraken_taxonomy_setup.smk
#
# This file provides a set of modular rules to download all required taxonomy
# source files from NCBI and GTDB, and then merge them to create a unified
# taxdump (nodes.dmp, names.dmp, merged.dmp). This unified taxdump is essential
# for the dual-database Kraken2 workflow.

# --- Configuration Helper ---
# Use a variable for the base path to keep rule definitions clean.
# This path is defined in the config.yaml file.
TAXONOMY_PATH = config["kraken_dual_db"]["taxonomy"]["path"]
GTDB_RELEASE = str(config["kraken_dual_db"]["taxonomy"]["gtdb_release"])

rule download_ncbi_taxonomy:
    """
    Downloads the NCBI taxonomy dump file (taxdmp.zip).
    """
    output:
        # We list the extracted files we need to ensure they are present
        nodes=os.path.join(TAXONOMY_PATH, "ncbi", "nodes.dmp"),
        names=os.path.join(TAXONOMY_PATH, "ncbi", "names.dmp"),
        merged=os.path.join(TAXONOMY_PATH, "ncbi", "merged.dmp")
    params:
        ncbi_tax_url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",
        target_dir=os.path.join(TAXONOMY_PATH, "ncbi")
    log:
        "logs/download_ncbi_taxonomy.log"
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml") # Assumes wget and unzip are in this environment
    threads:
        1
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    shadow: "copy-minimal"
    message:
        "Downloading NCBI taxonomy source files."
    shell:
        """
        # Create the target directory
        mkdir -p {params.target_dir}

        # Download the zip file
        wget -q -c -O {params.target_dir}/taxdmp.zip {params.ncbi_tax_url} > {log} 2>&1

        # Extract only the required dmp files from the zip archive
        unzip -q -o -d {params.target_dir} {params.target_dir}/taxdmp.zip nodes.dmp names.dmp merged.dmp >> {log} 2>&1
        """

rule download_gtdb_taxonomy:
    """
    Downloads the GTDB taxonomy files for a specific release.
    """
    output:
        ar_tax=os.path.join(TAXONOMY_PATH, "gtdb", f"ar53_taxonomy_r{GTDB_RELEASE}.tsv.gz"),
        bac_tax=os.path.join(TAXONOMY_PATH, "gtdb", f"bac120_taxonomy_r{GTDB_RELEASE}.tsv.gz")
    params:
        base_url=f"https://data.ace.uq.edu.au/public/gtdb/data/releases/release{GTDB_RELEASE}/{GTDB_RELEASE}.0",
        target_dir=os.path.join(TAXONOMY_PATH, "gtdb")
    log:
        "logs/download_gtdb_taxonomy.log"
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml") # Assumes wget is in this environment
    threads:
        1
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    shadow: "copy-minimal"
    message:
        "Downloading GTDB taxonomy source files for release {GTDB_RELEASE}."
    shell:
        """
        mkdir -p {params.target_dir}
        wget -q -c -P {params.target_dir} {params.base_url}/ar53_taxonomy_r{GTDB_RELEASE}.tsv.gz > {log} 2>&1
        wget -q -c -P {params.target_dir} {params.base_url}/bac120_taxonomy_r{GTDB_RELEASE}.tsv.gz >> {log} 2>&1
        """

rule merge_taxonomies:
    """
    Merges the downloaded NCBI and GTDB taxonomies into a unified taxdump.
    This rule calls a python script that accepts all file paths as arguments.
    """
    input:
        ncbi_nodes=os.path.join(TAXONOMY_PATH, "ncbi", "nodes.dmp"),
        ncbi_names=os.path.join(TAXONOMY_PATH, "ncbi", "names.dmp"),
        ncbi_merged=os.path.join(TAXONOMY_PATH, "ncbi", "merged.dmp"),
        gtdb_ar_tax=os.path.join(TAXONOMY_PATH, "gtdb", f"ar53_taxonomy_r{GTDB_RELEASE}.tsv.gz"),
        gtdb_bac_tax=os.path.join(TAXONOMY_PATH, "gtdb", f"bac120_taxonomy_r{GTDB_RELEASE}.tsv.gz")
    output:
        # The final merged taxdump files are the key output
        nodes=os.path.join(TAXONOMY_PATH, "merged", "nodes.dmp"),
        names=os.path.join(TAXONOMY_PATH, "merged", "names.dmp"),
        merged=os.path.join(TAXONOMY_PATH, "merged", "merged.dmp")
    params:
        merge_script=os.path.join(SRCDIR, "multitax_merge_ncbi_gtdb_taxonomy.py"),
        # The script needs a directory for its output files
        output_dir=os.path.join(TAXONOMY_PATH, "merged")
    log:
        "logs/merge_taxonomies.log"
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    threads:
        getThreads(12, bigmem=False, high_mem_cluster=False)
    resources:
        runtime = LONG_RUNTIME,
        mem_mb = HIGH_MEM_MB # Assumes HIGH_MEM_MB is defined in your config/functions
    shadow: "copy-minimal"
    message:
        "Merging NCBI and GTDB taxonomies into a unified taxdump."
    shell:
        """
        # The script will be updated to accept file paths directly.
        # This makes the rule clean and portable.
        python3 {params.merge_script} \\
            --ncbi-taxdump-dir $(dirname {input.ncbi_nodes}) \\
            --gtdb-files {input.gtdb_ar_tax} {input.gtdb_bac_tax} \\
            --output-dir {params.output_dir} \\
            --debug > {log} 2>&1
        """
