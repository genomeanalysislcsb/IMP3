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

# def get_max_cores_man(wildcards, attempt):
#     mem_per_core_attempt = get_min_mem_man(wildcards, attempt)
#     max_cores = config['cluster_specs']['main_cluster']['big_mem_cores_max']
#     total_mem = math.floor(config['cluster_specs']['main_cluster']['big_mem_total'] * 0.95)
#     attempt_supported_cores = math.floor(total_mem / mem_per_core_attempt)
#     return max(min(max_cores, attempt_supported_cores), 1)

# Mantis: protein annotation
# NOTE: assumming that mantis has been setup
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
        mantis_params="-gff -km -da heuristic",
        mem_avail_gb = math.floor(float(MAIN_BIGMEMTOTAL) * to_gb_dict[config['cluster_specs']['mem_unit']] * 0.75)
    threads:
        # get_max_cores_man
        getThreads(128, bigmem=True, high_mem_cluster=False)
    resources:
        runtime = LONG_RUNTIME,
        # mem_mb = get_min_mem_man
        mem_mb = MAIN_BIGMEMCORE_MB
    message:
        "Mantis: annotate {input}"
    conda:
        os.path.join(ENVDIR, "IMP_mantis.yaml")
    shell:
        """
        # Check installation
        echo "### Checking Mantis installation ... ###" >> {log}
        (mantis  check  --mantis_config {input.config} >> {log} 2>&1 && echo "Mantis is installed correctly." >> {log}) \
          || (echo "Mantis is not installed correctly." >> {log} && mantis setup --mantis_config {input.config} --cores {threads} >> {log} 2>&1)

        # Use 10% less cores/threads to try to avoiud OOM errors
        mantis_cores=$(echo "{threads} * 0.75" | bc | cut -d '.' -f1)

        # Run Mantis
        rm -rf Analysis/annotation/intermediary/mantis* >> {log} 2>&1
        mantis run -i {input.faa} \
                   -o {params.odir} \
                   -mc {input.config} \
                   {params.mantis_params} \
                   -c $mantis_cores \
                   -m {params.mem_avail_gb} >> {log} 2>&1

        # Compress all output files
        echo "Compressing Mantis output files" >> {log}
        for f in {params.odir}/*; do
            if [ -f $f ]; then
                gzip $f >> {log} 2>&1
            fi
        done
        echo "Done." >> {log}
        """

# Mantis: reformat consensus output for easier downstream processing
# NOTE: Cannot use Mantis' GFF file as it does not contain contigs as sequence-records.
#       The GFF file is created from the tabular output.
rule mantis_reformat_consensus:
    input:
        contab=rules.mantis_annot.output.contab
    output:
        contab="Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.tsv"  # Used to be temp()
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    message:
        "Mantis: reformat {input}"
    log:
        "logs/analysis/mantis_reformat.log"
    shell:
        # split by "|", merge annotations (replace "\t" by "|")
        # NOTE: using "|" as it should be very unlikely to have it in the annotation text
        """
        echo "Reformatting Mantis consensus output" >> {log}
        paste <(zcat {input.contab} 2>> {log} | cut -d '|' -f1 2>> {log}  | sed 's/\\t$//' 2>> {log}) \
              <(zcat {input.contab} 2>> {log} | cut -d '|' -f2- 2>> {log} | sed 's/^\\t//;s/\\t/|/g' 2>> {log}) > {output.contab} 2>> {log} 
        echo "Done." >> {log}
        """


checkpoint makegff_mantis:
    input:
        gff="Analysis/annotation/bakta/annotation.filt.gff", # TODO: use rules.annotate.output.XXX when the outputs are named in that rule
        mantis=rules.mantis_reformat_consensus.output.contab,
    output:
        gff="Analysis/annotation/intermediary/mantis/consensus_annotation.reformatted.gff"  # used to be temp()
    params:
        sep="|", # separator used in mantis_reformat_consensus
        prefix="Analysis/annotation/intermediary/mantis/consensus_annotation.mantis"
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        # Ensure at least 10gb of mem
        max(int(math.ceil(50 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    message:
        "Mantis: Create GFF from {input}"
    log:
        "logs/analysis/mantis_make_gff.log"
    script:
        os.path.join(SRCDIR, "Analysis", "mantis_gff.py")
