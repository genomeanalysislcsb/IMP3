import math

include:
    "mantis_config.smk"

# TODO: Check why 4000 mb/core with > 42 cores (168G total) or 8000 mb/core with > 25 (200G total) cores cant submit with slurm but 12000 mb/core + 18 cores (216G total) works 
# Maybe the problem is with the cluster configuration: for batch `DefMemPerCPU=1750 MaxMemPerCPU=1750`
# def get_min_mem_man(wildcards, attempt):
#     mem_per_core = float(config['cluster_specs']['main_cluster']['big_mem_per_core'])
#     mem_unit = config['cluster_specs']['mem_unit']
#     mem_avail_mb = math.floor(mem_per_core * to_mb_dict[mem_unit])
#     return attempt * max(mem_avail_mb, 4000)

def get_max_cores_man(wildcards, attempt):
    cores = max(getThreads(96, bigmem=True, high_mem_cluster=False), 1)
    return max(cores / attempt, 1)

# Mantis: protein annotation
# NOTE: assumming that mantis has been setup
def set_ma_runtime(wildcards, attempt):
    if attempt == 1:
        return SHORT_RUNTIME
    elif attempt == 2:
        return SHORT_RUNTIME
    else:
        return MEDIUM_RUNTIME

rule mantis_annot:
    input:
        faa="Analysis/annotation/bakta/annotation.faa",
        config=rules.mantis_config.output.config_file,
    output:
        outtab="Analysis/annotation/intermediary/mantis/output_annotation.tsv.gz",
        inttab="Analysis/annotation/intermediary/mantis/integrated_annotation.tsv.gz",
        contab="Analysis/annotation/intermediary/mantis/consensus_annotation.tsv.gz",
        intgff="Analysis/annotation/intermediary/mantis/integrated_annotation.gff.gz",
        congff="Analysis/annotation/intermediary/mantis/consensus_annotation.gff.gz",
        keggmodules="Analysis/annotation/intermediary/mantis/kegg_modules.tsv.gz",
        kosample="Analysis/annotation/intermediary/mantis/sample_kos.tsv.gz",
    log:
        "logs/analysis/mantis.log"
    params:
        odir=lambda wildcards, output: os.path.dirname(output.outtab),
        # mantis_params="-gff -km -da heuristic",
        mantis_params=config['mantis']['params'],
        mem_avail_gb = math.floor(float(MAIN_BIGMEMTOTAL) * to_gb_dict[config['cluster_specs']['mem_unit']] * 0.55),
        # mem_avail_gb = math.floor(float(HIGH_MEM_BIGMEMTOTAL) * to_gb_dict[config['cluster_specs']['mem_unit']] * 0.75),
        node_tmp_dir="/dev/shm",
        sub_threads = 12
    threads:
        get_max_cores_man
        # getThreads(128, bigmem=True, high_mem_cluster=False)
    retries:
        3
    resources:
        runtime = set_ma_runtime,
        mem_mb = MAIN_BIGMEMCORE_MB
    message:
        "Mantis: annotate {input}"
    conda:
        os.path.join(ENVDIR, "IMP_mantis.yaml")
    shadow: "copy-minimal"
    shell:
        """
        # First we copy the mantis config file and the fasta file to the node's tmp directory
        echo "### Copying Mantis config file and fasta file to node's tmp directory ... ###" >> {log}
        rsync -a {input.config} {params.node_tmp_dir}/mantis_config.yaml >> {log} 2>&1
        rsync -a {input.faa} {params.node_tmp_dir}/annotation.faa >> {log} 2>&1

        log_memory() {{
            used=$(free -m | awk 'NR==2 {{print $3}}')
            total=$(free -m | awk 'NR==2 {{print $2}}')
            gb=$(awk "BEGIN {{printf \\"%.2f\\", $used/1024}}")
            echo "Current memory usage: $gb GB ($used MB / $total MB)" >> {log}
        }}


        # Check installation
        echo "### Checking Mantis installation ... ###" >> {log}
        (mantis  check  --mantis_config {input.config} >> {log} 2>&1 && echo "Mantis is installed correctly." >> {log}) \
        || (echo "Mantis is not installed correctly." >> {log} && mantis setup --mantis_config {input.config} --cores {threads} >> {log} 2>&1)
        
        # Use N less cores to reduce memory pressure
        if [ {threads} -gt {params.sub_threads} ]; then
            mantis_cores=$(echo "{threads} - {params.sub_threads}" | bc)
        else
            mantis_cores={threads}
        fi

        # Run Mantis
        rm -rf {params.odir} >> {log} 2>&1
        
        log_memory

        mantis run -i {params.node_tmp_dir}/annotation.faa \
                   -o {params.node_tmp_dir}/mantis \
                   -mc {params.node_tmp_dir}/mantis_config.yaml \
                   {params.mantis_params} \
                   -c $mantis_cores \
                   -m {params.mem_avail_gb} >> {log} 2>&1
        
        log_memory

        # Compress all output files
        echo "Compressing Mantis output files" >> {log}
        for f in {params.node_tmp_dir}/mantis/*; do
            if [ -f $f ]; then
                gzip $f >> {log} 2>&1
            fi
        done

        # Ensure out dir is created
        mkdir -p {params.odir} >> {log} 2>&1

        # Copy output files back to the working directory
        echo "### Copying Mantis output files back to the working directory ... ###" >> {log}
        # We need to get the correct group in case we want to write e.g. a specific project directory
        TARGET_GROUP=$(stat -c '%G' "{params.odir}")
        rsync -a --chown=:$TARGET_GROUP {params.node_tmp_dir}/mantis/* {params.odir} >> {log} 2>&1

        log_memory

        echo "Done." >> {log}
        """

# Mantis: reformat consensus output for easier downstream processing
# NOTE: Cannot use Mantis' GFF file as it does not contain contigs as sequence-records.
#       The GFF file is created from the tabular output.
rule mantis_reformat_consensus:
    input:
        contab=rules.mantis_annot.output.contab
    output:
        contab="Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.tsv"  # Used to be te mp()
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    message:
        "Mantis: reformat {input}"
    log:
        "logs/analysis/mantis_reformat.log"
    shadow: "copy-minimal"
    shell:
        # split by "|", merge annotations (replace "\t" by "|")
        # NOTE: using "|" as it should be very unlikely to have it in the annotation text
        """
        echo "Reformatting Mantis consensus output" >> {log}
        paste <(zcat {input.contab} 2>> {log} | cut -d '|' -f1 2>> {log}  | sed 's/\\t$//' 2>> {log}) \
              <(zcat {input.contab} 2>> {log} | cut -d '|' -f2- 2>> {log} | sed 's/^\\t//;s/\\t/|/g' 2>> {log}) > {output.contab} 2>> {log} 
        echo "Done." >> {log}
        """


def set_mgm_runtime(wildcards, attempt):
    if attempt == 1:
        return VERY_SHORT_RUNTIME
    elif attempt == 2:
        return SHORTISH_RUNTIME
    else:
        return MEDIUM_RUNTIME

checkpoint makegff_mantis:
    input:
        gff="Analysis/annotation/bakta/annotation.filt.gff", # TODO: use rules.annotate.output.XXX when the outputs are named in that rule
        mantis=rules.mantis_reformat_consensus.output.contab,
    output:
        gff="Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.gff",
        ref_gffs_dir = directory("Analysis/annotation/intermediary/mantis/ref_gff")
    params:
        sep="|", # separator used in mantis_reformat_consensus
        prefix="Analysis/annotation/intermediary/mantis/ref_gff/consensus_annotation.mantis"
    resources:
        runtime = set_mgm_runtime,
        mem_mb = MAIN_MEMCORE_MB
    retries:
        3
    threads:
        # Ensure at least 10gb of mem
        max(int(math.ceil(50 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    message:
        "Mantis: Create GFF from {input}"
    log:
        "logs/analysis/mantis_make_gff.log"
    shadow: "copy-minimal"
    script:
        os.path.join(SRCDIR, "Analysis", "mantis_gff.py")
