rule get_bakta_db:
    output:
        bakta_db = os.path.join(DBPATH, 'bakta/db/bakta.db')
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        4
    conda:
        os.path.join(ENVDIR, "IMP_bakta.yaml")
    log:
        os.path.join(ENVDIR, "logs/bakta_downloadDB.log")
    params:
        db = os.path.join(DBPATH, 'bakta')
    message:
        "bakta: Downloading bakta database."
    shell:
        """
        bakta_db download --output {params.db} \
                          --type full > {log} 2>&1 \
         && touch {output} >> {log} 2>&1
        """
rule setup_bakta_functions_from_fork:
    output:
        bakta_path = directory(os.path.join(ROOTDIR, "submodules/bakta")), 
        ready_check = os.path.join(ROOTDIR, 'submodules/bakta.ready')
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "IMP_bakta.yaml")
    log:
        os.path.join(ENVDIR, "logs/bakta_setup_functions.log")
    params:
        fork_url = "https://github.com/ohickl/bakta.git",
        fork_branch = "dynamic_db_fetch"
    message:
        "bakta: Setting up functions from fork."
    shadow: "copy-minimal"
    shell:
        """
        # Clone git repo to submodule dir
        echo "Cloning bakta fork ..." >> {log}
        git clone -b {params.fork_branch} {params.fork_url} {output.bakta_path} >> {log} 2>&1

        # Find the bakta package directory in the conda environment
        BAKTA_CONDA_DIR=$(find ${{CONDA_PREFIX}} -type d -name "bakta" -path "*/site-packages/bakta" | head -n 1)

        if [ -z "$BAKTA_CONDA_DIR" ]; then
            echo "Error: Bakta directory not found in conda environment" >> {log}
            exit 1
        fi

        echo "Replacing functions in conda env ..." >> {log}
        # Replace files in the features directory
        for file in nc_rna.py nc_rna_region.py r_rna.py t_rna.py tm_rna.py s_orf.py cds.py crispr.py; do
            cp {output.bakta_path}/bakta/features/$file ${{BAKTA_CONDA_DIR}}/features/ >> {log} 2>&1
        done

        # Replace files in the main bakta directory
        for file in main.py config.py utils.py psc.py pscc.py ips.py ups.py; do
            cp {output.bakta_path}/bakta/bakta/$file ${{BAKTA_CONDA_DIR}}/ >> {log} 2>&1
        done

        # Increase diamond block size from 3 to 6
        echo "Increasing diamond block size ..." >> {log}
        find ${{BAKTA_CONDA_DIR}} -type f -name '*.py' -exec sed -i "s/'--block-size', '3'/'--block-size', '6'/g" {{}} \; >> {log} 2>&1

        # Touch ready file
        touch {output.ready_check} >> {log} 2>&1

        echo "Done." >> {log}
        """

def set_ba_runtime(wildcards, attempt):
    if attempt == 1:
        return SHORTISH_RUNTIME
    elif attempt == 2:
        return SHORT_RUNTIME
    else:
        return MEDIUM_RUNTIME

rule annotate:
    input:
        db_ready = os.path.join(DBPATH, 'bakta/db/bakta.db'),
        bakta_check = os.path.join(ROOTDIR, 'submodules/bakta.ready'),
        assembly = 'Assembly/%s.assembly.merged.fa' % ASS
    output:
        bakta_dir = directory("Analysis/annotation/bakta"),
        gff_filt = "Analysis/annotation/bakta/annotation.filt.gff",
        faa = "Analysis/annotation/bakta/annotation.faa",
        fna = "Analysis/annotation/bakta/annotation.fna",
        ffn = "Analysis/annotation/bakta/annotation.ffn",
        gff_arch = "Analysis/annotation/bakta/annotation.gff3.gz"
    threads:
        getThreads(128, bigmem=False, high_mem_cluster=False)
    resources:
        runtime = set_ba_runtime,
        mem_mb = MAIN_MEMCORE_MB
    retries:
        3
    log:
        "logs/analysis/annotate.log"
    params:
        db = os.path.join(DBPATH, "bakta/db"),
        tmp = TMPDIR,
        prefix = "annotation",
        gff = "Analysis/annotation/bakta/annotation.gff3"
    conda:
        os.path.join(ENVDIR, "IMP_bakta.yaml")
    message:
        "annotate: Running bakta."
    shadow: "copy-minimal"
    shell:
        """
        echo "Running bakta ..." >> {log}
        bakta --db {params.db} \
              --tmp_db_path /dev/shm \
              --verbose \
              --debug \
              --prefix {params.prefix} \
              --keep-contig-headers \
              --meta \
              --force \
              --skip-write-genbank-embl \
              --tmp-dir {params.tmp} \
              --output {output.bakta_dir} \
              --threads {threads} \
              {input.assembly} >> {log} 2>&1
        
        # Filter gff to just contain annotation lines
        echo "Filtering gff file ..." >> {log}

        LN=`grep -Hn "^>" {params.gff} | head -n1 | cut -f2 -d ":" || if [[ $? -eq 141 ]]; then true; else exit $?; fi`
        LN1=1
        LN=$(($LN-$LN1))
        head -n $LN {params.gff} | grep -P -v "^#|\tBakta\tregion\t" | sort | uniq | grep -v "^==" > {output.gff_filt}

        # Compress gff
        echo "Compressing gff file ..." >> {log}
        pigz -p {threads} -c {params.gff} > {output.gff_arch} 2>> {log} \
         && rm -rf {params.gff} >> {log} 2>&1

        echo "Done." >> {log}
        """

