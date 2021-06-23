localrules: VIS

stats_figs =[]
if "preprocessing" in IMP_STEPS or MG or MT:
    stats_figs.append(expand("Visualization/{types}.reads.png",types=TYPES))
if "assembly" in IMP_STEPS or CONTIGS:
    stats_figs.append(expand("Visualization/Assembly_{fig}.png", fig=["length_hist","cov_hist","gc_hist","cov_gc_scatter"]))
if "analysis" in IMP_STEPS:
    stats_figs.append(expand("Visualization/Analysis_{fig}.png", fig=["length_hist","cov_hist","cov_scatter","func_sets"]))
if "taxonomy" in IMP_STEPS:
    stats_figs.append(expand("Visualization/Taxonomy_{type}_kraken_tree.png", type=TYPES))
    stats_figs.append("Visualization/Taxonomy_kraken_plot.png")
    stats_figs.append(expand("Visualization/Taxonomy_{type}_mOTU_tree.png", type=TYPES))
    stats_figs.append("Visualization/Taxonomy_mOTU_plot.png")
if "binning" in IMP_STEPS:    
    stats_figs.append("Visualization/Binning_alluvial.png")
    stats_figs.append(expand("Visualization/Binning_{type}_depth.png",type=TYPES))
    if "binny" in BIN_STEPS:
        stats_figs.append("Visualization/Binning_vizbin.png")
    
rule vis_stats:
    input:
        stats_figs
    output:
        touch("status/stats.visualized")

rule visualize_read_stats:
    input:
        "Stats/all_stats.Rdata"
    output:
        report(expand("Visualization/{types}.reads.png",types=TYPES),
               caption="../../report/stats/Reads.rst",category="Preprocessing")
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    log: "logs/visualize_read_stats.log"
    message: "visualize_read_stats: Plotting reads figure."
    threads: 1
    conda: ENVDIR + "/IMP_visualize.yaml"
    script:
        SRCDIR + "/visualize_reads.R"

rule visualize_assembly_stats:
    input:
        "Stats/all_stats.Rdata"
    output:
        report("Visualization/Assembly_length_hist.png",caption="../../report/stats/Assembly_length.rst",category="Assembly"),
        report("Visualization/Assembly_cov_hist.png",caption="../../report/stats/Assembly_cov.rst",category="Assembly"),
        report("Visualization/Assembly_gc_hist.png",caption="../../report/stats/Assembly_gc.rst",category="Assembly"),
        report("Visualization/Assembly_cov_gc_scatter.png",caption="../../report/stats/Assembly_cov_gc.rst",category="Assembly")
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    log: "logs/visualize_assembly_stats.log"
    message: "visualize_assembly_stats: Plotting figures for contigs."
    threads: 1
    conda: ENVDIR + "/IMP_visualize.yaml"
    script:
        SRCDIR + "/visualize_assembly.R"

rule visualize_analysis_stats:
    input:
        "Stats/all_stats.Rdata"
    output:
        report("Visualization/Analysis_length_hist.png",caption="../../report/stats/Analysis_length.rst",category="Analysis"),
        report("Visualization/Analysis_cov_hist.png",caption="../../report/stats/Analysis_cov.rst",category="Analysis"),
        report("Visualization/Analysis_cov_scatter.png",caption="../../report/stats/Analysis_cov_scat.rst",category="Analysis"),
        report("Visualization/Analysis_func_sets.png",caption="../../report/stats/Analysis_func.rst",category="Analysis")
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    log: "logs/visualize_analysis_stats.log"
    message: "visualize_analysis_stats: Plotting gene figures."
    threads: 1
    conda: ENVDIR + "/IMP_visualize.yaml"
    script:
        SRCDIR + "/visualize_analysis.R"
 
if CONTIGS or "assembly" in IMP_STEPS:        
    rule visualize_kraken_reads:
        input:
            "Stats/all_stats.Rdata",
            "Analysis/taxonomy/kraken/%s.kraken.parsed.tsv" %ASS
        output:
            report(expand("Visualization/Taxonomy_{type}_kraken_tree.png",type=TYPES),caption="../../report/stats/Taxonomy_kraken_withAssembly.rst",category="Taxonomy"),
            report("Visualization/Taxonomy_kraken_plot.png",caption="../../report/stats/Taxonomy_kraken_scatter_withAssembly.rst",category="Taxonomy")
        resources:
            runtime = "2:00:00",
            mem = MEMCORE
        log: "logs/visualize_kraken_reads.log"
        message: "visualize_kraken_reads: Plotting taxonomy results."
        threads: 1
        conda: ENVDIR + "/IMP_visualize.yaml"
        script:
            SRCDIR + "/visualize_taxonomy_kraken_withAssembly.R"
else:
    rule visualize_kraken_reads:
        input:
            "Stats/all_stats.Rdata"
        output:
            report(expand("Visualization/Taxonomy_{type}_kraken_tree.png",type=TYPES),caption="../../report/stats/Taxonomy_kraken.rst",category="Taxonomy"),
            report("Visualization/Taxonomy_kraken_plot.png",caption="../../report/stats/Taxonomy_kraken_scatter.rst",category="Taxonomy")
        resources:
            runtime = "2:00:00",
            mem = MEMCORE
        log: "logs/visualize_kraken_reads.log"
        message: "visualize_kraken_reads: Plotting taxonomy results."
        threads: 1
        conda: ENVDIR + "/IMP_visualize.yaml"
        script:
            SRCDIR + "/visualize_taxonomy_kraken.R"

if "analysis" in IMP_STEPS:
    rule visualize_mOTU_reads:
        input:
            "Stats/all_stats.Rdata",
            "Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv"
        output:
            report(expand("Visualization/Taxonomy_{type}_mOTU_tree.png",type=TYPES),caption="../../report/stats/Taxonomy_mOTU_withAnalysis.rst",category="Taxonomy"),
            report("Visualization/Taxonomy_mOTU_plot.png",caption="../../report/stats/Taxonomy_mOTU_withAnalysis.rst",category="Taxonomy"),
        resources:
            runtime = "2:00:00",
            mem = MEMCORE
        log: "logs/visualize_mOTU_reads.log"
        message: "visualize_mOTU_reads: Plotting taxonomy results."
        threads: 1
        conda: ENVDIR + "/IMP_visualize.yaml"
        script:
            SRCDIR + "/visualize_mOTU_withAnalysis.R"
else:
    rule visualize_mOTU_reads:
        input:
            "Stats/all_stats.Rdata"
        output:
            report(expand("Visualization/Taxonomy_{type}_mOTU_tree.png",type=TYPES),caption="../../report/stats/Taxonomy_mOTU.rst",category="Taxonomy"),
            report("Visualization/Taxonomy_mOTU_plot.png",caption="../../report/stats/Taxonomy_mOTU.rst",category="Taxonomy"),
        resources:
            runtime = "2:00:00",
            mem = MEMCORE
        log: "logs/visualize_mOTU_reads.log"
        message: "visualize_mOTU_reads: Plotting taxonomy results."
        threads: 1
        conda: ENVDIR + "/IMP_visualize.yaml"
        script:
            SRCDIR + "/visualize_mOTU.R"

if "binny" in BIN_STEPS:
    rule visualize_bins:
        input:
            stats="Stats/all_stats.Rdata",
            bins=expand("Binning/{binner}/scaffold2bin.tsv",binner=BIN_STEPS),
            dastool="Binning/selected_DASTool_scaffolds2bin.txt",
            vizbin="Binning/binny/%s.vizbin.with-contig-names.points" %ASS
        output:
            report("Visualization/Binning_alluvial.png",caption="../../report/stats/Binning_alluvial.rst",category="Binning"),
            report(expand("Visualization/Binning_{type}_depth.png",type=TYPES),caption="../../report/stats/Binning_depth.rst",category="Binning"),
            report("Visualization/Binning_vizbin.png",caption="../../report/stats/Binning_vizbin.rst",category="Binning")
        resources:
            runtime = "2:00:00",
            mem = MEMCORE
        log: "logs/visualize_bins.log"
        message: "visualize_bins: Plotting binning results."
        threads: 1
        conda: ENVDIR + "/IMP_visualize.yaml"
        script:
            SRCDIR + "/visualize_binning.R"
else:
    rule visualize_bins:
        input:
            stats="Stats/all_stats.Rdata",
            bins=expand("Binning/{binner}/scaffold2bin.tsv",binner=BIN_STEPS),
            dastool="Binning/selected_DASTool_scaffolds2bin.txt"
        output:
            report("Visualization/Binning_alluvial.png",caption="../../report/stats/Binning_alluvial.rst",category="Binning"),
            report(expand("Visualization/Binning_{type}_depth.png",type=TYPES),caption="../../report/stats/Binning_depth.rst",category="Binning")
        resources:
            runtime = "2:00:00",
            mem = MEMCORE
        log: "logs/visualize_bins.log"
        message: "visualize_bins: Plotting binning results."
        threads: 1
        conda: ENVDIR + "/IMP_visualize.yaml"
        script:
            SRCDIR + "/visualize_binning.R"




# master command
rule VIS:
    input:
        vis_target
    output:
        touch('status/visualize.done')
