db_output = ["DB","DB/contigs.bson"]
if "binning" in IMP_STEPS:
    db_output.append("DB/bins.bson")

localrules: prep_mongo, shut_mongo, mongo_server

rule mongo_server:
    input:
        db_output
    output:
        touch("status/DB_filling.running")
    threads: 1
    params:
        port = config["mongo_port"]
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    log: "logs/summary_mongo_server.log"
    message: "mongoDB_server: Running mongo."
    threads: 1
    conda: ENVDIR + "/IMP_multi.yaml"
    shell:
        """
        numactl --interleave=all mongod --dbpath {input[0]} --port {params.port} --auth --logpath {log}
        """

rule shut_mongo:
    input:
        db_output
    output:
        touch("status/DB_filling.done")
    params:
        port = config["mongo_port"],
    message: "shut_mongo: killing mongod on config["mongo_port"]."
    shell:
        """
        mongod --dbpath {input[0]} --port {params.port} --shutdown
        """

rule prep_mongo:
    input:
        "status/workflow.done"
    output:
        directory("DB")
    message: "prep_mongo: making DB directoy."
    shell:
        """
        mkdir -p {output}
        """

rule mongoDB_bins:
    input:
        "DB"
    output:
        pipe("DB/bins.bson")
    threads: 1
    params:
        ass = ASS,
        sample = SAMPLE,
        binners = config["Binning"]["Binners"],
        dbs = config["hmm_DBs"],
        user = config["mongo_user"],
        password = config["mongo_password"],
        port = config["mongo_port"]
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    log: "logs/summary_mongoDB_bins.log"
    message: "mongoDB_bins: Filling bin documents."
    threads: 1
    conda: ENVDIR + "/IMP_multi.yaml"
    script:
        SRCDIR + "/mongify_bins_single_IMP3.py"
         
rule mongoDB_contigs:
    input:
        "DB"
    output:
        pipe("DB/contigs.bson")
    threads: 1
    params:
        ass = ASS,
        sample = SAMPLE,
        binners = config["Binning"]["Binners"],
        dbs = config["hmm_DBs"],
        user = config["mongo_user"],
        password = config["mongo_password"],
        port = config["mongo_port"],
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    log: "logs/summary_mongoDB_contigs.log"
    message: "mongoDB_bins: Filling contig documents."
    threads: 1
    conda: ENVDIR + "/IMP_multi.yaml"
    script:
        SRCDIR + "/mongify_contigs_single_IMP3.py"


        

