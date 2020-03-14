PLOT_OUTPUT_FILES = [
    "Analysis/results/IMP-vizbin_length.png",
    "Analysis/results/IMP-vizbin_length_GC.png",
    "Analysis/results/IMP-vizbin_standard.png",
    "Analysis/results/IMP-mg_reads_density.png",
    "Analysis/results/IMP-mg_coverage_density.png",
    "Analysis/results/IMP-mg_vizbin_length_cov.png",
    "Analysis/results/IMP-mg_vizbin_length_depth.png",
    "Analysis/results/IMP-mg_var_count.png",
    "Analysis/results/IMP-mg_vizbin_length_vardens.png",
    "Analysis/results/mg.read_stats.html",
    "Analysis/results/mg.read_stats.txt",
    "Analysis/results/mg_mapping_stats.html",
    "Analysis/results/mg_mapping_stats.txt",
    "Analysis/results/mg_results.Rdat",
    expand("Binning/binny/reachabilityDistanceEstimates.{pk}.{nn}.tsv \
    Binning/binny/clusterFirstScan.{pk}.{nn}.tsv \
    Binning/binny/bimodalClusterCutoffs.{pk}.{nn}.tsv \
    Binning/binny/contigs2clusters.{pk}.{nn}.tsv \
    Binning/binny/contigs2clusters.{pk}.{nn}.RDS \
    Binning/binny/clusteringWS.{pk}.{nn}.Rdata \
    #AHB: Binning/binny/finalClusterMap.{pk}.{nn}.pdf".split(),pk=config["binning"]["pk"],nn=config["binning"]["nn"])
    Binning/binny/finalClusterMap.{pk}.{nn}.pdf".split(),pk=config["binning"]["pk"],nn=config["binning"]["nn"])
]



rule visualize:
    input:
       'Analysis/mg.read_counts.txt',
       'Analysis/mg.assembly.contig_flagstat.txt',
       'Analysis/mg.assembly.contig_coverage.txt',
       'Analysis/mg.assembly.contig_depth.txt',
       'Analysis/mg.variants.samtools.vcf.gz',
       'Analysis/mg.assembly.gc_content.txt',
       "Analysis/mg.vizbin.with-contig-names.points",
       "Analysis/annotation/annotation_CDS_RNA_hmms.gff" #to be adjusted for gaia
#       "Analysis/results/quast/summary/combined_quast_output/contigs_reports/nucmer_output/aux/mgmt.assembly.merged.coords_edited"
    output:
        PLOT_OUTPUT_FILES
    shell:
        """
        PLOT_SCRIPT="{SRCDIR}/IMP_visualize_mg_commented.R"

        echo "[x] PLOT `date +"%Y/%m/%d %H:%M:%S"`"
        mkdir -p Analysis/results
        mkdir -p Binning/binny/clusterFiles
	#AHB: Rscript $PLOT_SCRIPT Analysis/results {input} {SRCDIR}/IMP_plot_binny_functions.R {config[binning][pk]} {config[binning][nn]}
        Rscript $PLOT_SCRIPT Analysis/results {input} {SRCDIR}/IMP_plot_binny_functions.R {config[binning][pk]} {config[binning][nn]}
        """
# config is adjusted to gaia style
