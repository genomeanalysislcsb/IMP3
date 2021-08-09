rule GRID_summarize:
    input:
        grid=gather_multi_bins,
        dastool=expand("Dereplication/{sample}/selected_DASTool_summary.txt",
                 sample=samples.loc[[a and b for a,b in zip(samples.ass == 1,samples.rebinning == "yes")],"sample"])
    output:
        "Dereplication/dRep_multi/allSamples_bin_stats.tsv"
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        sample = samples.loc[[a and b for a,b in zip(samples.ass == 1,samples.rebinning == "yes")],"sample"].tolist()
    log: "logs/dereplication_GRID_summarize.log"
    message: "GRID_summarize: Summarizing binning results."
    threads: 1
    conda: ENVDIR + "/IMP_binning.yaml"
    script:
        SRCDIR + "/multi_collect_GRID.R"
    

rule getMultiBins:
    input: 
        bins="Dereplication/{sample}/selected_DASTool_summary.txt",
        indir="Dereplication/{sample}/selected_DASTool_bins"
    output:
        "status/dRep_multi_{sample}.dRep_multi.prepped"
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    message: "getMultiBins: Copying DAStool bins of {wildcards.sample} to common directory and to single directories."
    run:
        selected_bins=pd.read_csv(input.bins,sep="\t")
        for BIN in selected_bins['bin']:
            shell("cp {input.indir}/{BIN}.contigs.fa Dereplication/dRep_multi/{wildcards.sample}.{BIN}.contigs.fa")
            shell("mkdir -p Dereplication/GRiD_multi/{wildcards.sample}.{BIN} && cp Dereplication/{wildcards.sample}/selected_DASTool_bins/{BIN}.contigs.fa Dereplication/GRiD_multi/{wildcards.sample}.{BIN}/contigs.fa")
        shell("touch {output}")


localrules: pre_grid
checkpoint pre_grid:
    input:
       expand("status/dRep_multi_{sample}.dRep_multi.prepped",
              sample=samples.loc[[a and b for a,b in zip(samples.ass == 1,samples.rebinning == "yes")],"sample"])
    output:
       touch("status/dereplication_GRiD_multi_grid.prepped")

rule grid:
    input:
       "Dereplication/GRiD_multi/{clusterID}/contigs.fa",
       "Dereplication/GRiD_multi/{clusterID}/reads.fastq.gz"
    output:
        "Dereplication/GRiD_multi/{clusterID}/grid/reads.GRiD.txt",
        "Dereplication/GRiD_multi/{clusterID}/grid/reads.pdf"
    params:
        clusterInput= "Dereplication/GRiD_multi/{clusterID}",
        outdir = "Dereplication/GRiD_multi/{clusterID}/grid"
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_grid.yaml"
    log: "logs/dereplication_grid.{clusterID}.log"
    message: "grid: Running GRiD on {wildcards.clusterID}."
    shell:
        """
        grid single -r {params.clusterInput} -e fastq.gz -o {params.outdir} -g {input[0]} > {log} 2>&1
        """

rule prepare_gridBed:
    input:
        "Dereplication/GRiD_multi/{clusterID}/contigs.fa"
    output:
        temporary("Dereplication/GRiD_multi/{clusterID}/contigs.bed")
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    threads: 1
    log: "logs/dereplication_prepare_gridBed.{clusterID}.log"
    message: "prepare_gridBed: Preparing .bed file for GRiD - {wildcards.clusterID}."
    shell:
        """
        {SRCDIR}/fastaBed.pl {input} >> {output} 2> {log}
        """

rule prepare_gridReads:
    input:
        lambda wildcards: expand("{raw_directory}/Assembly/mg.reads.sorted.bam",
                raw_directory=samples.path[wildcards.sample]),
        "Dereplication/GRiD_multi/{sample}.{clusterID}/contigs.bed"
    output:
        temporary("Dereplication/GRiD_multi/{sample}.{clusterID}/reads.fastq.gz")
    wildcard_constraints:
        sample="[^.]*"
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: 1
    log: "logs/dereplication_prepare_gridReads.{sample}.{clusterID}.log"
    message: "prepare_gridReads: Preparing fastq file for GRiD - {wildcards.sample}{wildcards.clusterID}."
    conda: ENVDIR + "/IMP_mapping.yaml"
    shell:
        """
        samtools view -b -F 5 -L {input[1]} {input[0]} 2>> {log} | samtools fastq -n - > {output} 2>> {log} 
        samtools view -b -f 65 -F 4 -L {input[1]} {input[0]} 2>> {log} | samtools fastq -n - >> {output} 2>> {log}
        samtools view -b -f 129 -F 4 -L {input[1]} {input[0]} 2>> {log} | samtools fastq -n - >> {output} 2>> {log}
        """
