localrules: dRep_simple_dir

rule dRep_simple_dir:
    input:
        "sample_table.tsv"
    output:
        directory("Dereplication/dRep_single")
    shell:
        """
        mkdir -p {output}
        """

rule prep_dRep_single:
    input:
        dir="Dereplication/dRep_single",
        bins= lambda wildcards: expand("{raw_directory}/Binning/selected_DASTool_summary.txt",
                raw_directory=samples.path[wildcards.sample]),
        indir= lambda wildcards: expand("{raw_directory}/Binning/selected_DASTool_bins",
                raw_directory=samples.path[wildcards.sample])
    output:
        "status/dereplication_{sample}.dRep_single.prepped"
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    threads: 1
    message: "prep_dRep_single: Copying DAStool bins of {wildcards.sample} to common directory."
    run:
        selected_bins=pd.read_csv(input.bins[0],sep="\t")
        for BIN in selected_bins['bin']:
            shell("cp {input.indir}/{BIN}.contigs.fa Dereplication/dRep_single/{wildcards.sample}.{BIN}.contigs.fa")
        shell("touch {output}")

rule dRep_simple:
    input:
        ctrl=expand("status/dereplication_{sample}.dRep_single.prepped",
         sample=samples.loc[samples.bins == 1,"sample"]),
        dir="Dereplication/dRep_single"
    output:
        directory("Dereplication/dRep_single_out/dereplicated_genomes"),
        expand("Dereplication/dRep_single_out/data_tables/{file}.csv",
                file=["Chdb","Cdb","Sdb","Wdb","Widb","Ndb"])
    threads: getThreads(BIGCORENO)
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    params:
        out = "Dereplication/dRep_single_out",
        db = config['db_path'],
        length = config['dereplication']['dRep']['length'],
        completeness = config['dereplication']['dRep']['completeness'],
        contamination = config['dereplication']['dRep']['contamination'],
        checkM_method = config['dereplication']['dRep']['checkM_method'],
        MASH_sketch = config['dereplication']['dRep']['MASH_sketch'],
        S_algorithm = config['dereplication']['dRep']['S_algorithm'],
        P_ani = config['dereplication']['dRep']['P_ani'],
        S_ani = config['dereplication']['dRep']['S_ani'],
        SkipMash = config['dereplication']['dRep']['SkipMash'],
        SkipSecondary = config['dereplication']['dRep']['SkipSecondary'],
        cov_thresh = config['dereplication']['dRep']['cov_thresh'],
        coverage_method = config['dereplication']['dRep']['coverage_method'],
        clusterAlg = config['dereplication']['dRep']['clusterAlg'],
        completeness_weight = config['dereplication']['dRep']['completeness_weight'],
        contamination_weight = config['dereplication']['dRep']['contamination_weight'],
        strain_heterogeneity_weight = config['dereplication']['dRep']['strain_heterogeneity_weight'],
        N50_weight = config['dereplication']['dRep']['N50_weight'],
        size_weight = config['dereplication']['dRep']['size_weight']
    log: "logs/dereplication_dRep_simple.log"
    message: "dRep_simple: Running dRep on bins of input samples."
    conda: ENVDIR + "/IMP_multi.yaml"
    shell:
        """
        checkm data setRoot {params.db}/checkM/
        dRep filter {params.out} -d -p {threads} --length {params.length} --completeness {params.completeness} \
         --contamination {params.contamination} --checkM_method {params.checkM_method} \
         -g {input.dir}/*contigs.fa >> {log} 2>&1
        dRep cluster {params.out} -p {threads} --MASH_sketch {params.MASH_sketch} --S_algorithm {params.S_algorithm} \
         --P_ani {params.P_ani} --S_ani {params.S_ani} {params.SkipMash} \
         {params.SkipSecondary} --cov_thresh {params.cov_thresh} \
         --coverage_method {params.coverage_method} --clusterAlg {params.clusterAlg} >> {log} 2>&1
        dRep choose {params.out} -p {threads} --completeness_weight {params.completeness_weight} \
         --contamination_weight {params.contamination_weight} --strain_heterogeneity_weight {params.strain_heterogeneity_weight} \
         --N50_weight {params.N50_weight} --size_weight {params.size_weight} --checkM_method {params.checkM_method} >> {log} 2>&1
        dRep evaluate {params.out} --evaluate all -p {threads} >> {log} 2>&1
        """
        
