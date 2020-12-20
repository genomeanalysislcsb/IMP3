localrules: prepare_binny

rule prepare_binny:
    input:
       mgdepth='Stats/mg/%s.assembly.contig_depth.txt' % ASS,
       vizbin='Binning/binny/%s.vizbin.with-contig-names.points' % ASS,
       gff='Analysis/annotation/annotation_CDS_RNA_hmms.gff'
    output:
       directory("Binning/binny/clusterFiles")
#    threads: 1
    message: "Prepare binny."
    shell:
       """
       mkdir -p {output} || echo "{output} exists"
       """

rule binny:
    input:
       outdir="Binning/binny/clusterFiles",
       mgdepth='Stats/mg/%s.assembly.contig_depth.txt' % ASS,
       vizbin="Binning/binny/%s.vizbin.with-contig-names.points" % ASS,
       gff="Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        expand("Binning/binny/reachabilityDistanceEstimates.{pk}.{nn}.tsv \
        Binning/binny/clusterFirstScan.{pk}.{nn}.tsv \
        Binning/binny/bimodalClusterCutoffs.{pk}.{nn}.tsv \
        Binning/binny/contigs2clusters.{pk}.{nn}.tsv \
        Binning/binny/contigs2clusters.{pk}.{nn}.RDS \
        Binning/binny/clusteringWS.{pk}.{nn}.Rdata".split(),pk=config["binning"]["binny"]["pk"],nn=config["binning"]["binny"]["nn"]),
        report(expand("Binning/binny/finalClusterMap.{pk}.{nn}.png",pk=config["binning"]["binny"]["pk"],nn=config["binning"]["binny"]["nn"]),caption="../../../report/stats/binny_final_map.rst",category="Binning"),
        "Binning/binny/binny_WS.Rdata"
    params:
        plot_functions = SRCDIR + "/IMP_plot_binny_functions.R",
        binnydir="Binning/binny"
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_binning.yaml" 
    log: "logs/binning_binny.log"
    message: "binny: Running Binny."
    script:
        SRCDIR + "/binny.R"

rule tar_binny_files:
    input:
        "Binning/binny/clusterFiles",
        expand("Binning/binny/contigs2clusters.{pk}.{nn}.tsv",pk=config["binning"]["binny"]["pk"],nn=config["binning"]["binny"]["nn"])
    output:
        "Binning/binny/clusterFiles.tar.gz"
    threads: 1
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    params:
        intermediary = "Binning/binny/clusterFiles/"
    log: "logs/binning_tar_binny_files.log"
    message: "tar_binny_files: Compressing intermediary files from binny."
    shell:
       """
       tar cvzf {output} {params.intermediary} >> {log} 2>&1 && rm -r {params.intermediary} >> {log} 2>&1
       """
