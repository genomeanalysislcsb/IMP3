rule GTDBtk_bins:
    input:
        "Dereplication/dRep_multi_out/dereplicated_genomes/{clusterID}.contigs.fa"
    output:
        directory("Dereplication/GTDB_multi/{clusterID}.GTDB_out")
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    params:
        outdir = "Dereplication/GTDB_multi/{clusterID}"
    conda: ENVDIR + "/IMP_gtdbtk.yaml"
    log: "logs/dereplication_GTDBtk_bins.{clusterID}.log"
    message: "GTDBtk_bins: Classifiying {wildcards.clusterID} using GTDBtk."
    threads: getThreads(BIGCORENO)
    shell:
        """
        mkdir -p {params.outdir} && cp {input} {params.outdir}
        export GTDBTK_DATA_PATH="{DBPATH}/GTDB_tk"
        export PYTHONPATH=$CONDA_PREFIX/lib/python3.7/site-packages
        gtdbtk classify_wf --genome_dir {params.outdir}/ -x fa --out_dir {output[0]} --cpus {threads} > {log} 2>&1
        """

rule GTDBtk_summarize:
    input:
        gtdb = gather_multi_GTDB_bins,
        dRep = expand("Dereplication/dRep_multi_out/data_tables/{file}.csv",
                file=["Chdb","Cdb","Sdb","Wdb","Widb","Ndb"]),
        bins= "Dereplication/dRep_multi/allSamples_bin_stats.tsv"
    output:
        "Dereplication/dRep_multi_out/dereplicated_rebin_stats.tsv",
        report("Dereplication/rebinned_bin_stats.tsv",category="Re-binning and dereplication")
    resources:
        runtime = "4:00:00",
        mem = MEMCORE
    log: "logs/dereplication_GTDBtk_summarize.log"
    message: "GTDBtk_summarize: Summarizing binning results."
    threads: 1
    conda: ENVDIR + "/IMP_binning.yaml"
    script:
        SRCDIR + "/multi_summarize_GTDBtk.R"

