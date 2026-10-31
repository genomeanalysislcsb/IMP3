import math

localrules: ALL_stats

gathered = [report("Stats/all_stats.tsv",category="Stats"),
            "Stats/all_stats.Rdata"]
            
if "binning" in IMP_STEPS:
    gathered.append(report("Stats/all_bin_stats.tsv",category="Binning"))

import os
import math

bracken_levels_cfg = config['bracken']['levels'].split()

def check_step(step):
    return step in IMP_STEPS

# Define helper functions for conditional inputs
def get_annotation_files(wildcards):
    ann = config["annotation"]
    files = []
    if ann == "mantis":
        checkpoint_output = checkpoints.makegff_mantis.get(**wildcards).output[0]
        refs = glob_wildcards('Analysis/annotation/intermediary/mantis/ref_gff/consensus_annotation.mantis.{i}.gff').i
        for t in TYPES:
            files += expand("Analysis/annotation/{type}.mantis.{ref}.counts.tsv", type=TYPES, ref=refs)
            files += expand("Analysis/annotation/{type}.mantis.{ref}.counts.tsv.summary", type=TYPES, ref=refs)
    elif ann == "hmmer":
        hmmer_dbs = config["hmm_DBs"].split()
        for t in TYPES:
            for db in hmmer_dbs:
                files.append(f"Analysis/annotation/{t}.hmmer.{db}.counts.tsv.summary")
    return files

def get_binning_contigs(wildcards):
    checkpoint_output = checkpoints.getBins.get(**wildcards).output[0]
    return expand("Binning/MetaWrap/bins/{bin}/contigs.bed", bin=glob_wildcards('Binning/MetaWrap/bins/{i}/contigs.fasta').i) + ["Binning/MetaWrap/bins.contigs"]

rule gather_stats:
    input:
        # Control Files
        wf_control="status/workflow.done",
        prep_control=lambda wildcards: "status/preprocessing.analysed" if not check_step("preprocessing") else [],
        cont_control=lambda wildcards: "status/assembly.analysed" if CONTIGS else [],
        
        # Statistics and Reports
        stats_files = expand("Stats/{type}/{type}.read_stats.tsv", type=TYPES),
        krak_reps = expand("Analysis/taxonomy/kraken/{type}.kraken.report", type=TYPES),
        bracken = expand("Analysis/taxonomy/bracken/{type}/{level}.bracken", type=TYPES, level=bracken_levels_cfg) if check_step("taxonomy") and MG else [],
        motus_counts = expand("Analysis/taxonomy/mOTUs/{type}.mOTU.counts.tsv", type=TYPES) if check_step("taxonomy") else [],
        
        # Annotation Files
        annotation_files = get_annotation_files,

        # Assembly Stats
        assembly_contig_flagstat = expand("Stats/{type}/{ass}.assembly.contig_flagstat.txt", type=TYPES, ass=ASS),
        assembly_length = expand("Stats/{ass}/{ass}.assembly.length.txt", ass=ASS),
        assembly_gc_content = expand("Stats/{ass}/{ass}.assembly.gc_content.txt", ass=ASS),
        assembly_contig_depth = expand("Stats/{type}/{ass}.assembly.contig_depth.txt", type=TYPES, ass=ASS),
        assembly_breadth_thresholds = expand("Stats/{type}/{ass}.assembly.thresholds.bed.gz", type=TYPES, ass=ASS),
        
        # Kraken Assembly Parsed
        kraken_assembly_parsed = expand("Analysis/taxonomy/kraken/{ass}.kraken.contig.unfiltered.report", ass=ASS),
        
        # Annotation GFF Files
        gene_gff_mantis = "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        gene_gff_filt = "Analysis/annotation/bakta/annotation.filt.gff",
        
        # Gene Length and Depth
        gene_length = expand("Stats/{type}/{type}.gene.len", type=TYPES),
        gene_depth = expand("Stats/{type}/{type}.gene_depth.avg", type=TYPES),
        
        # Binning Stats
        binning_stats = expand("Stats/{ass}/{ass}.bins.tsv", ass=ASS) if 'binning' in IMP_STEPS else [],
        binning_contigs = get_binning_contigs if 'binning' in IMP_STEPS else [],
        
        # mOTU Links
        mOTU_links = "Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv"
    output:
        gathered = "Stats/all_stats.tsv",
        gathered_rdat = "Stats/all_stats.Rdata",
        gathered_bins = "Stats/all_bin_stats.tsv" if 'binning' in IMP_STEPS else [],
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        steps = IMP_STEPS,
        ass = ASS,
        types = TYPES,
        ann = config["annotation"],
        hmmer_dbs = config['hmm_DBs'].split() if config['annotation'] == "hmmer" else "",
        mantis_dir = "Analysis/annotation/" if config["annotation"] == "mantis" else ""
    log:
        "logs/misc/gather_stats.log"
    message:
        "gather_stats: Summarizing stats."
    threads:
        # Ensure at least 5gb of mem
        max(int(math.ceil(5 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    shadow: "copy-minimal"
    script:
        SRCDIR + "/summarize_all_stats.R"
         
rule ALL_stats:
    input:
        gathered
    output:
        touch("status/stats.done")

