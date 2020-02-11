PLOT_OUTPUT_FILES = [
    "Stats/mg/mg.read_stats.html",
    "Stats/mg/mg.read_stats.txt",
    "Stats/mg/mg_mapping_stats.html",
    "Stats/mg/mg_mapping_stats.txt",
    "Stats/mt/mt.read_stats.html",
    "Stats/mt/mt.read_stats.txt",
    "Stats/mt/mt_mapping_stats.html",
    "Stats/mt/mt_mapping_stats.txt",
    "Stats/mgmt/mgmt_results.Rdat"
]

rule visualize_mgmt:
    input:
       'Analysis/mg.read_counts.txt',
       'Analysis/mt.read_counts.txt',
       'Analysis/mg.assembly.contig_flagstat.txt',
       'Analysis/mt.assembly.contig_flagstat.txt',
       'Analysis/mg.assembly.contig_coverage.txt',
       'Analysis/mt.assembly.contig_coverage.txt',
       'Analysis/mg.assembly.contig_depth.txt',
       'Analysis/mt.assembly.contig_depth.txt',
       'Analysis/mgmt.assembly.gc_content.txt',
       "Analysis/mgmt.vizbin.with-contig-names.points",
       "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        PLOT_OUTPUT_FILES
    shell:
        """
        PLOT_SCRIPT="{SRCDIR}/IMP_visualize_mgmt_commented.R"

        echo "[x] PLOT `date +"%Y/%m/%d %H:%M:%S"`"
        mkdir -p Analysis/results
        mkdir -p Binning/binny/clusterFiles
        Rscript $PLOT_SCRIPT {config[imp_lib]} Analysis/results {input} {SRCDIR}/IMP_plot_binny_functions.R {config[binning][binny][pk]} {config[binning][binny][nn]}
        """


