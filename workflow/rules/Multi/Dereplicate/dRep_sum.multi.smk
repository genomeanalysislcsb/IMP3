rule dRep_sum_multi:
    input:
        dRep = expand("Dereplication/dRep_multi_out/data_tables/{file}.csv",
                file=[Chdb","Cdb","Sdb","Wdb","Widb","Ndb"]),
        bins= "Dereplication/dRep_multi/allSamples_bin_stats.tsv"
    output:
        "Dereplication/dRep_multi_out/dereplicated_rebin_stats.tsv",
        report("Dereplication/rebinned_bin_stats.tsv",category="Re-binning and dereplication")
    threads: 1 
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_multi.yaml"
    log: "logs/dereplication_dRep_sum_multi.log"
    message: "dRep_sum_multi: Gathering data on dereplicated bins."
    script:
        SRCDIR + "/multi_rebin_derepOutput.R" 
