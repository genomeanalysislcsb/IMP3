rule dRep_sum_simple:
    input:
        dRep = expand("Dereplication/dRep_single_out/data_tables/{file}.csv",
                file=["Chdb","Cdb","Sdb","Wdb","Widb","Ndb"]),
        ori = expand("{raw_directory}/Stats/{ass}/{ass}.bins.tsv",
                raw_directory=samples.loc[samples.bins == 1,"path"],
                ass=samples.loc[samples.bins == 1,"assembly_type"]) 
    output:
        "Dereplication/dRep_single_out/originalNdereplicated_bin_stats.tsv",
        report("Dereplication/dereplicated_bin_stats.tsv",category="Re-binning and dereplication")
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        sample = samples.loc[samples.bins == 1,"sample"].tolist()
    conda: ENVDIR + "/IMP_multi.yaml"
    log: "logs/dereplication_dRep_sum_simple.log"
    message: "dRep_sum_simple: Gathering data on dereplicated bins."
    script:
        SRCDIR + "/multi_derepOutput.R" 
