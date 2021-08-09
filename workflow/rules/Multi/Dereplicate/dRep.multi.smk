localrules: dRep_dir

rule dRep_dir:
    input:
        "sample_table.tsv"
    output:
        directory("Dereplication/dRep_multi")
    shell:
        """
        mkdir -p {output}
        """

checkpoint dRep_multi:
    input:
        ctrl=expand("status/dRep_multi_{sample}.dRep_multi.prepped",
          sample=samples.loc[[a and b for a,b in zip(samples.ass == 1,samples.rebinning == "yes")],"sample"]),
    output:
        directory("Dereplication/dRep_multi_out/dereplicated_genomes"),
        expand("Dereplication/dRep_multi_out/data_tables/{file}.csv",
                file=["Chdb","Cdb","Sdb","Wdb","Widb","Ndb"])
    threads: getThreads(BIGCORENO)
    conda: ENVDIR + "/IMP_multi.yaml"
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    params:
        dir = "Dereplication/dRep_multi",
        out = "Dereplication/dRep_multi_out",
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
    log: "logs/dereplication_dRep_multi.log"
    message: "dRep_multi: Running dRep on multi-sample bins."
    shell:
        """
        checkm data setRoot {params.db}/checkM/
        dRep filter {params.out} -p {threads} --length {params.length} --completeness {params.completeness} \
          --contamination {params.contamination} --checkM_method {params.checkM_method} \
          -g {params.dir}/*contigs.fa >> {log} 2>&1
        dRep cluster {params.out} -p {threads} --MASH_sketch {params.MASH_sketch} --S_algorithm {params.S_algorithm} \
          --P_ani {params.P_ani} --S_ani {params.S_ani} {params.SkipMash}  \
          {params.SkipSecondary} --cov_thresh {params.cov_thresh} \
          --coverage_method {params.coverage_method} --clusterAlg {params.clusterAlg} >> {log} 2>&1
        dRep choose {params.out} -p {threads} --completeness_weight {params.completeness_weight} \
          --contamination_weight {params.contamination_weight} --strain_heterogeneity_weight {params.strain_heterogeneity_weight} \
          --N50_weight {params.N50_weight} --size_weight {params.size_weight} --checkM_method {params.checkM_method} >> {log} 2>&1
        dRep evaluate {params.out} -e a -p {threads} >> {log} 2>&1
        """

