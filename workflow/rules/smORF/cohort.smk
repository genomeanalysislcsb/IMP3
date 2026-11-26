# All samples together (the paper's smORF_prediction_with_clusters): smORFs clustered with MMseqs2 at each level of
# clustering: levels (identity, coverage), representative genes, presence per sample and group, length
# distributions, numbers per sample and per level

def _metadata_groups():
    """sample -> group from metadata (first two columns; header line allowed); 'all' without metadata"""
    groups = {s: "all" for s in SAMPLES}
    md = config.get("metadata") or ""
    if md:
        with open(os.path.expandvars(md)) as fh:
            for line in fh:
                p = re.split(r"[\t,]", line.rstrip("\n"))
                if len(p) >= 2 and p[0] in groups and p[1].strip():
                    groups[p[0]] = p[1].strip()
    return groups


rule groups_csv:
    output:
        "Input/groups.csv"
    params:
        lines = "\n".join("%s,%s" % (s, g) for s, g in _metadata_groups().items())
    shell:
        "printf '%s\\n' '{params.lines}' > {output}"


rule concatenate_smorfs:
    input:
        expand("smORF/{s}.smorf.faa", s=SAMPLES)
    output:
        temp("TMP/all_smORF.fasta")
    threads: 1
    resources:
        runtime = 60,
        mem_mb = MEM
    shell:
        "cat {input} > {output}"


rule cluster_smorfs:
    input:
        "TMP/all_smORF.fasta"
    output:
        rep = "Cluster_predicted_smORF/Coverage_{cov}_Ident_{idt}_rep_seq.fasta",
        tsv = "Cluster_predicted_smORF/Coverage_{cov}_Ident_{idt}_cluster.tsv",
        all = "Cluster_predicted_smORF/Coverage_{cov}_Ident_{idt}_all_seqs.fasta"
    params:
        d = "Cluster_predicted_smORF/Coverage_{cov}_Ident_{idt}",
        tmp = os.path.join(TMP, "cluster_{cov}_{idt}"),
        p = MMSEQS_PARAMS
    threads: T_BIG
    resources:
        runtime = 720,
        mem_mb = BIGMEM * T_BIG
    log:
        "logs/Cluster_predicted_smORF/Coverage_{cov}_Ident_{idt}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "cluster_smorfs: smORFs of all samples, identity {wildcards.idt}, coverage {wildcards.cov}"
    shell:
        """
        mmseqs easy-cluster {input} {params.d} {params.tmp} -c {wildcards.cov} --min-seq-id {wildcards.idt} \\
            --remove-tmp-files 1 {params.p} --threads {threads} -v 1 > {log} 2>&1
        rm -rf {params.tmp}
        """


rule fromAA_to_DNA_predicted:
    input:
        g = expand("smORF/{s}.smorf.fna", s=SAMPLES),
        p = "Cluster_predicted_smORF/{lvl}_rep_seq.fasta"
    output:
        "Cluster_predicted_smORF/{lvl}_rep_seq.fna"
    params:
        script = os.path.join(SRCDIR, "fromAAtoDNAseqs_predicted.py")
    threads: T_SMALL
    resources:
        runtime = 300,
        mem_mb = MEM * T_SMALL
    log:
        "logs/fromAA_to_DNA_predicted/{lvl}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "fromAA_to_DNA_predicted: genes of the representative smORFs ({wildcards.lvl})"
    shell:
        "python {params.script} -g smORF -p {input.p} -o {output} -n {threads} -e .smorf.fna 2> {log}"


rule presence_absence_predicted_smORF:
    input:
        i = "Cluster_predicted_smORF/{lvl}_cluster.tsv",
        g = "Input/groups.csv"
    output:
        "Cluster_predicted_smORF/{lvl}_rep_presence_in_groups_report.tsv"
    params:
        script = os.path.join(SRCDIR, "smORF_presence_in_samples.py")
    threads: 4
    resources:
        runtime = 240,
        mem_mb = BIGMEM * 4
    log:
        "logs/presence_absence/{lvl}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "presence_absence_predicted_smORF: samples per group of each representative smORF ({wildcards.lvl})"
    shell:
        "python {params.script} -i {input.i} -g {input.g} -o {output} -t {threads} > {log} 2>&1"


rule predicted_size_distribution:
    input:
        s = expand("smORF/{s}.smorf.faa", s=SAMPLES),
        g = "Input/groups.csv"
    output:
        "Size_distribution/smORF/Per_Group_size_distribution.png"
    params:
        script = os.path.join(SRCDIR, "SEQ_len_dis_plot_julia.py")
    threads: 1
    resources:
        runtime = 60,
        mem_mb = BIGMEM
    log:
        "logs/size_distribution/per_group.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        "MPLBACKEND=Agg python {params.script} -i smORF --group_info {input.g} -o {output} > {log} 2>&1"


rule predicted_size_distribution_cluster:
    input:
        fasta = "Cluster_predicted_smORF/{lvl}_rep_seq.fasta",
        g = "Input/groups.csv"
    output:
        "Size_distribution/smORF/{lvl}_cluster_representative_size_distribution.png"
    params:
        script = os.path.join(SRCDIR, "SEQ_len_dis_plot_cluster_julia.py")
    threads: 1
    resources:
        runtime = 60,
        mem_mb = BIGMEM
    log:
        "logs/size_distribution/{lvl}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        "MPLBACKEND=Agg python {params.script} -i {input.fasta} --group_info {input.g} -o {output} > {log} 2>&1"


rule smorf_cohort_summary:
    input:
        samples = expand("Stats/{s}.smorf_stats.tsv", s=SAMPLES),
        tsv = expand("Cluster_predicted_smORF/{lvl}_cluster.tsv", lvl=LEVELS),
        rep = expand("Cluster_predicted_smORF/{lvl}_rep_seq.fasta", lvl=LEVELS)
    output:
        samples = "smORF_samples.tsv",
        clusters = "smORF_clusters.tsv",
        rep_samples = expand("Cluster_predicted_smORF/{lvl}_rep_samples.tsv", lvl=LEVELS)
    params:
        script = os.path.join(SRCDIR, "smorf_summary.py"),
        levels = " ".join("%s,Cluster_predicted_smORF/%s_cluster.tsv,Cluster_predicted_smORF/%s_rep_seq.fasta"
                          % (l, l, l) for l in LEVELS)
    threads: 1
    resources:
        runtime = 120,
        mem_mb = BIGMEM
    log:
        "logs/smorf_cohort_summary.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        """
        python {params.script} cohort --samples {input.samples} --clusters {params.levels} \\
            --out {output.samples} --cluster-summary {output.clusters} > {log} 2>&1
        """
