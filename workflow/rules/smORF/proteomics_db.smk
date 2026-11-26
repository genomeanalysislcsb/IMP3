# Per-sample metaproteomics search databases (the paper's MP_DB): proteins of the IMP3 run (Bakta) and the smORFs
# of the sample, the human RefSeq proteome and the cRAP contaminants, each step deduplicated with MMseqs2 (100 %
# identity and coverage); identical proteins of different samples get one ID (proteomics_db/protein_clusters.tsv)
#   proteomics_db/final/<sample>.faa      search database of the sample
#   proteomics_db/protein_clusters.tsv    cluster, membership, protein ID, source (imp / smorf / human / crap),
#                                         sample, contig, gene, smORF (<= 100 aa)

MPDB_PREFIX = MPDB.get("id_prefix", "imp3")
MPDB_HUMAN = _true(MPDB.get("human", True))
for s in SAMPLES:
    if not imp_file(s, "bakta_faa") or not imp_file(s, "bakta_tsv"):
        raise Exception("smORF module: proteomics_db needs the Bakta proteins and table (%s, %s) of %s"
                        % (INPUTS["bakta_faa"], INPUTS["bakta_tsv"], s))


rule mpdb_reference_headers:
    input:
        human = HUMAN_FAA if MPDB_HUMAN else [],
        crap = CRAP_FAA
    output:
        human = "proteomics_db/reference/human.faa",
        crap = "proteomics_db/reference/crap.faa"
    params:
        script = os.path.join(SRCDIR, "mpdb_headers.py"),
        prefix = MPDB_PREFIX
    threads: 1
    resources:
        runtime = 60,
        mem_mb = MEM * 4
    log:
        "logs/proteomics_db/reference_headers.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        """
        if [[ -n "{input.human}" ]]; then
            python {params.script} human --fasta {input.human} --prefix {params.prefix} --out {output.human} 2> {log}
        else
            : > {output.human}
        fi
        python {params.script} crap --fasta {input.crap} --prefix {params.prefix} --out {output.crap} 2>> {log}
        """


rule mpdb_human_crap:
    input:
        human = "proteomics_db/reference/human.faa",
        crap = "proteomics_db/reference/crap.faa"
    output:
        "proteomics_db/reference/human_crap_rep_seq.fasta"
    params:
        d = "proteomics_db/reference/human_crap",
        tmp = os.path.join(TMP, "mpdb_human_crap")
    threads: T_SMALL
    resources:
        runtime = 120,
        mem_mb = MEM * T_SMALL
    log:
        "logs/proteomics_db/human_crap.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "mpdb_human_crap: human proteome and cRAP without duplicates"
    shell:
        """
        mmseqs easy-cluster {input.human} {input.crap} {params.d} {params.tmp} -c 1.0 --min-seq-id 1.0 \\
            --cov-mode 0 --cluster-mode 0 --threads {threads} -v 1 > {log} 2>&1
        rm -rf {params.tmp}
        """


rule mpdb_sample_headers:
    input:
        faa = lambda wc: imp_file(wc.s, "bakta_faa"),
        tsv = lambda wc: imp_file(wc.s, "bakta_tsv"),
        smorf = "smORF/{s}.smorf.faa"
    output:
        imp = "proteomics_db/imp/{s}.faa",
        smorf = "proteomics_db/smorf/{s}.faa"
    params:
        script = os.path.join(SRCDIR, "mpdb_headers.py"),
        prefix = MPDB_PREFIX,
        max_aa = int(PRED.get("max_aa", 100))
    threads: 1
    resources:
        runtime = 60,
        mem_mb = MEM * 2
    log:
        "logs/proteomics_db/headers_{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        """
        python {params.script} imp --fasta {input.faa} --tsv {input.tsv} --sample {wildcards.s} \\
            --prefix {params.prefix} --max-len {params.max_aa} --out {output.imp} 2> {log}
        python {params.script} smorf --fasta {input.smorf} --sample {wildcards.s} --prefix {params.prefix} \\
            --out {output.smorf} 2>> {log}
        """


rule mpdb_sample:
    input:
        imp = "proteomics_db/imp/{s}.faa",
        smorf = "proteomics_db/smorf/{s}.faa",
        ref = "proteomics_db/reference/human_crap_rep_seq.fasta"
    output:
        "proteomics_db/per_sample/{s}.faa"
    params:
        d1 = "proteomics_db/imp_smorf/{s}",
        d2 = "proteomics_db/per_sample/{s}",
        tmp = os.path.join(TMP, "mpdb_{s}")
    threads: T_SMALL
    resources:
        runtime = 240,
        mem_mb = MEM * T_SMALL
    log:
        "logs/proteomics_db/sample_{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "mpdb_sample: IMP3 proteins + smORFs, then + human / cRAP, without duplicates ({wildcards.s})"
    shell:
        """
        mkdir -p $(dirname {params.d1})
        mmseqs easy-cluster {input.imp} {input.smorf} {params.d1} {params.tmp}_1 -c 1.0 --min-seq-id 1.0 \\
            --cov-mode 0 --cluster-mode 0 --threads {threads} -v 1 > {log} 2>&1
        mmseqs easy-cluster {params.d1}_rep_seq.fasta {input.ref} {params.d2} {params.tmp}_2 -c 1.0 --min-seq-id 1.0 \\
            --cov-mode 0 --cluster-mode 0 --threads {threads} -v 1 >> {log} 2>&1
        mv {params.d2}_rep_seq.fasta {output}
        rm -rf {params.tmp}_1 {params.tmp}_2 {params.d2}_all_seqs.fasta {params.d1}_all_seqs.fasta
        """


rule mpdb_all_samples:
    input:
        expand("proteomics_db/per_sample/{s}.faa", s=SAMPLES)
    output:
        tsv = "proteomics_db/all_samples_cluster.tsv"
    params:
        d = "proteomics_db/all_samples",
        cat = os.path.join(TMP, "mpdb_all_samples.faa"),
        tmp = os.path.join(TMP, "mpdb_all")
    threads: T_BIG
    resources:
        runtime = 1440,
        mem_mb = BIGMEM * T_BIG
    log:
        "logs/proteomics_db/all_samples.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "mpdb_all_samples: identical proteins across the samples"
    shell:
        """
        mkdir -p $(dirname {params.cat})
        cat {input} > {params.cat}
        mmseqs easy-cluster {params.cat} {params.d} {params.tmp} -c 1.0 --min-seq-id 1.0 --cov-mode 0 \\
            --cluster-mode 0 --threads {threads} -v 1 > {log} 2>&1
        rm -rf {params.tmp} {params.cat} {params.d}_all_seqs.fasta {params.d}_rep_seq.fasta
        """


rule mpdb_ids:
    input:
        tsv = "proteomics_db/all_samples_cluster.tsv",
        faa = expand("proteomics_db/per_sample/{s}.faa", s=SAMPLES)
    output:
        tsv = "proteomics_db/protein_clusters.tsv",
        faa = expand("proteomics_db/final/{s}.faa", s=SAMPLES)
    params:
        script = os.path.join(SRCDIR, "mpdb_cluster_ids.py")
    threads: min(T_SMALL, 8)
    resources:
        runtime = 240,
        mem_mb = BIGMEM * 4
    log:
        "logs/proteomics_db/ids.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "mpdb_ids: one ID per identical protein in the per-sample databases"
    shell:
        """
        python {params.script} -tsv {input.tsv} -f {input.faa} -o {output.tsv} -d proteomics_db/final \\
            -t {threads} > {log} 2>&1
        """
