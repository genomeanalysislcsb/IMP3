localrules: vis_stats

rule vis_stats:
    input:
        "Visualization/Binning_mg_depth.png",
        "Visualization/Binning_vizbin.png"
    output:
        touch("status/stats.visualized")

rule visualize_bins:
    input:
        stats="Stats/all_stats.Rdata",
        bins="Binning/binny/scaffold2bin.tsv",
        vizbin="Binning/binny/%s.vizbin.with-contig-names.points" %ASS
    output:
        report("Visualization/Binning_mg_depth.png",caption="../../report/stats/Binning_depth.rst",category="Binning"),
        report("Visualization/Binning_vizbin.png",caption="../../report/stats/Binning_vizbin.rst",category="Binning")
    params:
        runtime = "2:00:00",
        mem = MEMCORE
    log: "logs/visualize_bins.log"
    message: "visualize_bins: Plotting binning results."
    threads: 1
    conda: ENVDIR + "/IMP_visualize.yaml"
    script:
        SRCDIR + "/visualize_binny.R"


