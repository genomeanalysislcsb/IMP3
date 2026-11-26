# False positives among the smORF candidates (the paper's remove_false_positives.smk): ncRNA genes (Rfam, cmscan
# with the trusted cut-offs), tRNA / tmRNA genes (Aragorn, standard and bacterial code) and spurious ORFs (AntiFam,
# gathering cut-offs). A switched-off filter (prediction: rfam / antifam / aragorn: false) writes an empty result.

rule rfam_on_smorfs:
    input:
        db = expand(RFAM_CM + "{ext}", ext=[".i1f", ".i1i", ".i1m", ".i1p"]) if USE_RFAM else [],
        s = "potential_smORF/{s}.smorf.fna"
    output:
        "potential_smORF/RFAM/{s}.tbl"
    params:
        cm = RFAM_CM,
        use = "yes" if USE_RFAM else "no"
    threads: T_BIG
    resources:
        runtime = 600,
        mem_mb = MEM * T_BIG
    log:
        "logs/rfam_on_smorfs/{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "rfam_on_smorfs: Rfam on the smORF genes of {wildcards.s}"
    shell:
        """
        if [[ {params.use} == yes ]]; then
            cmscan --rfam --fmt 2 --nohmmonly --cut_tc -g --cpu {threads} --noali --tblout {output} -o /dev/null \\
                {params.cm} {input.s} 2> {log}
        else
            echo "rfam: false" > {log}; : > {output}
        fi
        """


rule antifam:
    input:
        s = "potential_smORF/{s}.smorf.faa",
        db = ANTIFAM_HMM if USE_ANTIFAM else []
    output:
        "potential_smORF/antifam/{s}.hits.tsv"
    params:
        script = os.path.join(SRCDIR, "antifam_search.py"),
        hmm = ANTIFAM_HMM,
        use = "yes" if USE_ANTIFAM else "no"
    threads: T_SMALL
    resources:
        runtime = 120,
        mem_mb = MEM * T_SMALL
    log:
        "logs/antifam/{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "antifam: spurious ORFs among the smORFs of {wildcards.s}"
    shell:
        """
        if [[ {params.use} == yes ]]; then
            python {params.script} -hmm {params.hmm} -seqs {input.s} -n {threads} -o {output} 2> {log}
            touch {output}
        else
            echo "antifam: false" > {log}; : > {output}
        fi
        """


rule aragorn:
    input:
        "potential_smORF/{s}.smorf.fna"
    output:
        "potential_smORF/Aragorn_{code}/{s}.fasta"
    wildcard_constraints:
        code = "std|gcbact"
    params:
        use = "yes" if USE_ARAGORN else "no"
    threads: 1
    resources:
        runtime = 960,
        mem_mb = MEM
    log:
        "logs/Aragorn_{code}/{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "aragorn: tRNA / tmRNA genes among the smORFs of {wildcards.s} (-{wildcards.code})"
    shell:
        """
        if [[ {params.use} == yes ]]; then
            aragorn -{wildcards.code} -d -w -fasta -o {output} {input} > {log} 2>&1
        else
            echo "aragorn: false" > {log}; : > {output}
        fi
        """


rule remove_false_positives:
    input:
        p = "potential_smORF/{s}.smorf.faa",
        g = "potential_smORF/{s}.smorf.fna",
        a = "potential_smORF/{s}.smorf.gff",
        rfam = "potential_smORF/RFAM/{s}.tbl",
        antifam = "potential_smORF/antifam/{s}.hits.tsv",
        arstd = "potential_smORF/Aragorn_std/{s}.fasta",
        arbac = "potential_smORF/Aragorn_gcbact/{s}.fasta"
    output:
        p = "smORF/{s}.smorf.faa",
        g = "smORF/{s}.smorf.fna",
        a = "smORF/{s}.smorf.gff"
    params:
        script = os.path.join(SRCDIR, "get_smORFs_sequences.py")
    threads: 1
    resources:
        runtime = 60,
        mem_mb = MEM * 2
    log:
        "logs/remove_false_positives/{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "remove_false_positives: final smORFs of {wildcards.s}"
    shell:
        """
        python {params.script} -p {input.p} -g {input.g} -a {input.a} -rfam {input.rfam} -antifam {input.antifam} \\
            -arstd {input.arstd} -arbac {input.arbac} -prot {output.p} -gene {output.g} -gff {output.a} > {log} 2>&1
        """


rule smorf_sample_stats:
    input:
        genes_default = "pyrodigal_gv/default/{s}.faa",
        genes_small = "pyrodigal_gv/small/{s}.faa",
        cand_default = "pyrodigal_gv/default/{s}.smorf.faa",
        cand_small = "pyrodigal_gv/small/{s}.smorf.faa",
        potential = "potential_smORF/{s}.smorf.faa",
        rfam = "potential_smORF/RFAM/{s}.tbl",
        antifam = "potential_smORF/antifam/{s}.hits.tsv",
        arstd = "potential_smORF/Aragorn_std/{s}.fasta",
        arbac = "potential_smORF/Aragorn_gcbact/{s}.fasta",
        final = "smORF/{s}.smorf.faa",
        imp3 = "pyrodigal_gv/imp3/{s}.smorf.faa" if INTERGENIC else [],
        regions = "Stats/{s}.intergenic.tsv" if INTERGENIC else []
    output:
        "Stats/{s}.smorf_stats.tsv"
    params:
        script = os.path.join(SRCDIR, "smorf_summary.py"),
        extra = lambda wc, input: ("--cand-imp3 %s --regions %s" % (input.imp3, input.regions)) if INTERGENIC else ""
    threads: 1
    resources:
        runtime = 30,
        mem_mb = MEM
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        """
        python {params.script} sample --sample {wildcards.s} --genes-default {input.genes_default} \\
            --genes-small {input.genes_small} --cand-default {input.cand_default} --cand-small {input.cand_small} \\
            --potential {input.potential} --rfam {input.rfam} --antifam {input.antifam} --aragorn-std {input.arstd} \\
            --aragorn-bact {input.arbac} --final {input.final} --out {output} {params.extra}
        """


rule smorf_metap_validation:
    input:
        smorfs = "smORF/{s}.smorf.faa",
        peptides = lambda wc: metap_file(wc.s) or [],
        imp3 = lambda wc: imp_file(wc.s, "bakta_faa") or [],
        human = HUMAN_FAA if (BUILD_MPDB and _true(MPDB.get("human", True))) or os.path.exists(HUMAN_FAA) else [],
        crap = CRAP_FAA
    output:
        tab = "validation/{s}.smorf_metap.tsv",
        stats = "validation/{s}.metap_stats.tsv"
    params:
        script = os.path.join(SRCDIR, "smorf_metap.py"),
        peptides = lambda wc, input: "--peptides %s" % input.peptides if input.peptides else "",
        max_aa = int(PRED.get("max_aa", 100)),
        min_pep = int(VAL.get("min_unique_peptides", 1))
    threads: 1
    resources:
        runtime = 120,
        mem_mb = MEM * 4
    log:
        "logs/validation/{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "smorf_metap_validation: smORFs of {wildcards.s} with identified peptides"
    shell:
        """
        python {params.script} --sample {wildcards.s} --smorfs {input.smorfs} {params.peptides} \\
            --background {input.imp3} {input.human} {input.crap} --max-aa {params.max_aa} \\
            --min-peptides {params.min_pep} --out {output.tab} --stats {output.stats} > {log} 2>&1
        """


rule annotation_with_smorfs:
    input:
        imp3 = lambda wc: HAVE_GFF[wc.s],
        gff = "smORF/{s}.smorf.gff",
        faa = "smORF/{s}.smorf.faa",
        status = "validation/{s}.smorf_metap.tsv",
        reads = lambda wc: "read_support/%s.smorf_reads.tsv" % wc.s if wc.s in RS_SAMPLES else []
    output:
        gff = "smORF/{s}.annotation_with_smORFs.gff",
        table = "smORF/{s}.smorf_table.tsv"
    params:
        script = os.path.join(SRCDIR, "smorf_gff.py"),
        reads = lambda wc, input: "--reads %s" % input.reads if input.reads else ""
    threads: 1
    resources:
        runtime = 60,
        mem_mb = MEM * 4
    log:
        "logs/annotation_with_smorfs/{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "annotation_with_smorfs: IMP3 annotation of {wildcards.s} with its smORFs (new file), table of its smORFs"
    shell:
        """
        python {params.script} --imp3-gff {input.imp3} --smorf-gff {input.gff} --smorf-faa {input.faa} \\
            --status {input.status} {params.reads} --sample {wildcards.s} --out {output.gff} \\
            --table {output.table} > {log} 2>&1
        """


rule smorf_tables_all:
    input:
        tables = expand("smORF/{s}.smorf_table.tsv", s=SAMPLES),
        stats = expand("validation/{s}.metap_stats.tsv", s=SAMPLES)
    output:
        table = "smORF_table.tsv",
        stats = "validation_samples.tsv"
    params:
        script = os.path.join(SRCDIR, "smorf_summary.py")
    threads: 1
    resources:
        runtime = 60,
        mem_mb = MEM
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        """
        python {params.script} concat --inputs {input.tables} --out {output.table}
        python {params.script} concat --inputs {input.stats} --out {output.stats}
        """
