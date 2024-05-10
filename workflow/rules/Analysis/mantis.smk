import math

include:
    "mantis_config.smk"

# TODO: Check why 4000 mb/core + 56 cores or 8000 mb/core + 28 cores cores cant submit with slurm but manually. 12000 + 16 works 

def get_min_mem_man(wildcards, attempt):
    mem_per_core = float(config['cluster_specs']['main_cluster']['big_mem_per_core'])
    mem_unit = config['cluster_specs']['mem_unit']
    mem_avail_mb = math.floor(mem_per_core * to_mb_dict[mem_unit])
    return attempt * max(mem_avail_mb, 4000)

def get_max_cores_man(wildcards, attempt):
    mem_per_core_attempt = get_min_mem_man(wildcards, attempt)
    max_cores = config['cluster_specs']['main_cluster']['big_mem_cores_max']
    total_mem = math.floor(config['cluster_specs']['main_cluster']['big_mem_total'] * 0.95)
    attempt_supported_cores = math.floor(total_mem / mem_per_core_attempt)
    return max(min(max_cores, attempt_supported_cores), 1)

# Mantis: protein annotation
# NOTE: assumming that mantis has been setup
rule mantis_annot:
    input:
        faa="Analysis/annotation/prokka.faa",
        config=rules.mantis_config.output.config_file,
    output:
        outtab="Analysis/annotation/intermediary/mantis/output_annotation.tsv",
        inttab="Analysis/annotation/intermediary/mantis/integrated_annotation.tsv",
        contab="Analysis/annotation/intermediary/mantis/consensus_annotation.tsv",
        intgff="Analysis/annotation/intermediary/mantis/integrated_annotation.gff",
        congff="Analysis/annotation/intermediary/mantis/consensus_annotation.gff",
        keggmodules="Analysis/annotation/intermediary/mantis/kegg_modules.tsv",
        kosample="Analysis/annotation/intermediary/mantis/sample_kos.tsv",
    log:
        "logs/analysis/mantis.log"
    params:
        odir=lambda wildcards, output: os.path.dirname(output.outtab),
        params="-gff -km -da heuristic",
    threads:
        get_max_cores_man
    resources:
        runtime = LONG_RUNTIME,
        mem_mb = get_min_mem_man
    retries:
        3
    message:
        "Mantis: annotate {input}"
    conda:
        os.path.join(ENVDIR, "IMP_mantis.yaml")
    shell:
        """
        rm -rf Analysis/annotation/intermediary/mantis*
        mantis run -i {input.faa} \
                   -o {params.odir} \
                   -mc {input.config} \
                   {params.params} \
                   -c {threads} &> {log}
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
    shell:
        # split by "|", merge annotations (replace "\t" by "|")
        # NOTE: using "|" as it should be very unlikely to have it in the annotation text
        "paste <(cut -d '|' -f1 {input.contab} | sed 's/\\t$//') <(cut -d '|' -f2- {input.contab} | sed 's/^\\t//;s/\\t/|/g') > {output.contab}"

checkpoint makegff_mantis:
    input:
        gff="Analysis/annotation/annotation.filt.gff", # TODO: use rules.annotate.output.XXX when the outputs are named in that rule
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
