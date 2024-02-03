# SortMeRNA:
#   https://github.com/biocore/sortmerna/
#   https://github.com/biocore/sortmerna/wiki/User-manual-v4.0

# NOTE
#   The helper function below should be moved later to a designated file.
#
#   The check whether the working directory exists and is not empty
#   should not be required (as snakemake should remove it when (re-)running the rule)
#   but it is included anyway for additional safety.

def smk_get_sortmerna_params(wildcards, input):
    """For snakemake: get SortMeRNA parameters"""
    # reference files
    # --fastx: output aligned reads
    # -m: memory requirements for indexing
    # --aligned/--other: output (non)-aligned reads
    params = "{refs} --fastx -m {mem} --aligned --other".format(
        refs=" ".join(["--ref %s" % ref for ref in input.refs]),
        mem=int(MAIN_MEMCORE[:-1]) * 1000
    )
    # additional parameters for PE reads
    if hasattr(input, "r1"):
        assert hasattr(input, "r2"), f"Found R1 but not R2 in {input}"
        # --otu-map: OTU map (input to QIIME's make_otu_table.py)
        # --paired_in: flags the PE reads as "aligned", when either of them is aligned
        # --out2: output paired reads into separate files
        params = f"{params} --otu_map --paired_in --out2"
    return params

rule filter_rrna_pe:
    input:
        r1="Preprocessing/mt.r1.trimmed.fq",
        r2="Preprocessing/mt.r2.trimmed.fq",
        refs=expand(
            "{dbpath}/sortmerna/{fasta}.fasta",
            dbpath=DBPATH, fasta=config["sortmerna"]["files"]
        ),
    output:
        filt_r1="Preprocessing/mt.r1.trimmed__rna_filtered.fq",
        filt_r2="Preprocessing/mt.r2.trimmed__rna_filtered.fq",
        rrna_r1="Preprocessing/mt.r1.trimmed.rna.fq",
        rrna_r2="Preprocessing/mt.r2.trimmed.rna.fq",
        otu="Preprocessing/otu_map.txt",
        # temp: working and index directories
        wdir=temp(directory("Preprocessing/sortmerna_pe")),
        idx=temp(directory("Preprocessing/sortmerna_idx")), # store index files separately so they can be re-used in filter_rrna_se
    log:
        "logs/preprocessing_filter_rrna_pe.log"
    threads:
        # Ensure at least 100gb of mem
        max(int(math.ceil(100 / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    params:
        params=smk_get_sortmerna_params,
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "filter_rna: filtering rRNA from mt PE reads."
    shell:
        """
        if [ -d "{output.wdir}" ] && [ "$(ls -A {output.wdir})" ]; then
            echo "ERROR: {output.wdir} exists and is not empty. Please remove it and re-run the pipeline." >> {log} 2>&1
            exit 1
        fi
        sortmerna {params.params} --threads {threads} --reads {input.r1} --reads {input.r2} --workdir {output.wdir}/ --idx-dir {output.idx}/ >> {log} 2>&1
        cat {output.wdir}/out/aligned.log >> {log}
        mv {output.wdir}/out/other_fwd.fq {output.filt_r1} >> {log} 2>&1
        mv {output.wdir}/out/other_rev.fq {output.filt_r2} >> {log} 2>&1
        mv {output.wdir}/out/aligned_fwd.fq {output.rrna_r1} >> {log} 2>&1
        mv {output.wdir}/out/aligned_rev.fq {output.rrna_r2} >> {log} 2>&1
        if [ -f {output.wdir}/out/otu_map.txt ]; then
            mv {output.wdir}/out/otu_map.txt {output.otu} >> {log} 2>&1
        else
            touch {output.otu}
            echo "# WARNING: SortMeRNA OTUs file does not exist" > status/sortmerna.warning
        fi
        """

rule filter_rrna_se:
    input:
        se="Preprocessing/mt.se.trimmed.fq",
        refs=expand(
            "{dbpath}/sortmerna/{fasta}.fasta",
            dbpath=DBPATH, fasta=config["sortmerna"]["files"]
        ),
        idx=rules.filter_rrna_pe.output.idx, # to avoid additional overhead
    output:
        filt_se="Preprocessing/mt.se.trimmed__rna_filtered.fq",
        rrna_se="Preprocessing/mt.se.trimmed.rna.fq",
        # temp: working directory
        wdir=temp(directory("Preprocessing/sortmerna_se")),
    log:
        "logs/preprocessing_filter_rrna_se.log"
    threads:
        # Ensure at least 100gb of mem
        max(int(math.ceil(100 / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    params:
        params=smk_get_sortmerna_params,
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "filter_rna: filtering rRNA from mt SE reads."
    shell:
        """
        if [ -d "{output.wdir}" ] && [ "$(ls -A {output.wdir})" ]; then
            echo "ERROR: {output.wdir} exists and is not empty. Please remove it and re-run the pipeline." >> {log} 2>&1
            exit 1
        fi
        sortmerna {params.params} --threads {threads} --reads {input.se} --workdir {output.wdir}/ --idx-dir {input.idx}/ >> {log} 2>&1
        cat {output.wdir}/out/aligned.log >> {log}
        mv {output.wdir}/out/other.fq {output.filt_se} >> {log} 2>&1
        mv {output.wdir}/out/aligned.fq {output.rrna_se} >> {log} 2>&1
        """

# If KEEP_RRNAs is set to True, then the rRNA-filtered reads are concatenated with rRNA reads
if KEEP_RRNAs:
    rule copy_sortmerna_out:
        input:
            r1_filt = rules.filter_rrna_pe.output.filt_r1,
            r2_filt = rules.filter_rrna_pe.output.filt_r2,
            se_filt = rules.filter_rrna_se.output.filt_se,
            r1_rrna = rules.filter_rrna_pe.output.rrna_r1,
            r2_rrna = rules.filter_rrna_pe.output.rrna_r2,
            se_rrna = rules.filter_rrna_se.output.rrna_se
        output:
            r1="Preprocessing/mt.r1.trimmed__sortmerna.fq",
            r2="Preprocessing/mt.r2.trimmed__sortmerna.fq",
            se="Preprocessing/mt.se.trimmed__sortmerna.fq"
        log:
            "logs/preprocessing_copy_sortmerna_out.log"
        threads:
            getThreads(4)
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        message:
            "filter_rna: concatenating rRNA reads with other reads."
        shell:
            """
            echo "Concatenating rRNA-filtered reads with rRNA reads ..." > {log}
            cat {input.r1_filt} {input.r1_rrna} > {output.r1} 2>> {log}
            cat {input.r2_filt} {input.r2_rrna} > {output.r2} 2>> {log}
            cat {input.se_filt} {input.se_rrna} > {output.se} 2>> {log}
            echo "Done." >> {log}
            """
else:
    rule copy_sortmerna_out:
        input:
            r1=rules.filter_rrna_pe.output.filt_r1,
            r2=rules.filter_rrna_pe.output.filt_r2,
            se=rules.filter_rrna_se.output.filt_se
        output:
            r1="Preprocessing/mt.r1.trimmed__sortmerna.fq",
            r2="Preprocessing/mt.r2.trimmed__sortmerna.fq",
            se="Preprocessing/mt.se.trimmed__sortmerna.fq"
        log:
            "logs/preprocessing_copy_sortmerna_out.log"
        threads:
            getThreads(4)
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        message:
            "filter_rna: copying rRNA-filtered reads."
        shell:
            """
            echo "Only copying rRNA-filtered reads ..." > {log}
            cp {input.r1} {output.r1} 2>> {log}
            cp {input.r2} {output.r2} 2>> {log}
            cp {input.se} {output.se} 2>> {log}
            echo "Done." >> {log}
            """
