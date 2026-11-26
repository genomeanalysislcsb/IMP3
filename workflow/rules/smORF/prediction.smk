# smORF prediction on the IMP3 assembly of each sample (rules of the paper's smORF_prediction.smk):
# two pyrodigal-gv runs (default: genes >= 60 nt; small: genes >= 15 nt), smORF candidates of each (complete genes
# with an RBS motif, <= max_aa), deduplicated with MMseqs2 (100 % identity and coverage)

rule link_contigs:
    input:
        lambda wc: imp_file(wc.s, "assembly")
    output:
        "Input/contigs/{s}.fasta.gz"
    shell:
        """
        if [[ $(head -c2 {input} | od -An -tx1 | tr -d ' ') == "1f8b" ]]; then
            ln -sf {input} {output}
        else
            gzip -c {input} > {output}
        fi
        """


rule pyrodigal:
    input:
        "Input/contigs/{s}.fasta.gz"
    output:
        a = "pyrodigal_gv/{run}/{s}.gff",
        d = "pyrodigal_gv/{run}/{s}.fna",
        p = "pyrodigal_gv/{run}/{s}.faa"
    wildcard_constraints:
        run = "default|small"
    params:
        mode = PRED.get("mode", "meta"),
        opts = lambda wc: PRED.get(wc.run, "--min-gene 60 --max-overlap 60" if wc.run == "default"
                                   else "--min-gene 15 --max-overlap 15")
    threads: T_SMALL
    resources:
        runtime = 240,
        mem_mb = MEM * T_SMALL
    log:
        "logs/pyrodigal_gv/{run}/{s}.log"
    benchmark:
        "benchmarks/pyrodigal_gv/{run}/{s}.tsv"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "pyrodigal: gene calling ({wildcards.run}) on {wildcards.s}"
    shell:
        """
        pyrodigal-gv -i {input} -p {params.mode} -d {output.d} -a {output.p} -f gff -o {output.a} -j {threads} \\
            -c -m --no-stop-codon {params.opts} > {log} 2>&1
        # gene, protein and contig names get the prefix of the run (default_ / small_)
        sed -i s/">"/">{wildcards.run}_"/g {output.p}
        sed -i s/">"/">{wildcards.run}_"/g {output.d}
        sed -i "/^[^#]/s/^/{wildcards.run}_/" {output.a}
        """


rule smorf_candidates:
    input:
        "pyrodigal_gv/{run}/{s}.faa"
    output:
        "pyrodigal_gv/{run}/{s}.smorf.faa"
    wildcard_constraints:
        run = "default|small"
    params:
        script = os.path.join(SRCDIR, "process_smorfs.py"),
        max_aa = int(PRED.get("max_aa", 100))
    threads: 1
    resources:
        runtime = 60,
        mem_mb = MEM
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "smorf_candidates: smORFs of {wildcards.s} ({wildcards.run})"
    shell:
        "python {params.script} -i {input} -o {output} -s {wildcards.s} -m {params.max_aa}"


rule deduplicate_smorfs:
    input:
        default = "pyrodigal_gv/default/{s}.smorf.faa",
        small = "pyrodigal_gv/small/{s}.smorf.faa"
    output:
        p = "potential_smORF/{s}.smorf.faa",
        tsv = "pyrodigal_gv/deduplicated_smorfs/{s}_cluster.tsv"
    params:
        d = "pyrodigal_gv/deduplicated_smorfs/{s}",
        tmp = os.path.join(TMP, "dedup_{s}")
    threads: T_BIG
    resources:
        runtime = 240,
        mem_mb = MEM * T_BIG
    log:
        "logs/deduplicate_smorfs/{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "deduplicate_smorfs: identical smORFs of the two runs on {wildcards.s} (MMseqs2)"
    shell:
        """
        mmseqs easy-cluster {input.default} {input.small} {params.d} {params.tmp} -c 1.0 --min-seq-id 1.0 \\
            --cov-mode 0 --cluster-mode 0 --threads {threads} -v 1 > {log} 2>&1
        mv {params.d}_rep_seq.fasta {output.p}
        rm -rf {params.tmp}
        """


rule get_DNA_sequences:
    input:
        p = "potential_smORF/{s}.smorf.faa",
        gd = "pyrodigal_gv/default/{s}.fna",
        gs = "pyrodigal_gv/small/{s}.fna",
        ad = "pyrodigal_gv/default/{s}.gff",
        asm = "pyrodigal_gv/small/{s}.gff"
    output:
        g = "potential_smORF/{s}.smorf.fna",
        a = "potential_smORF/{s}.smorf.gff"
    params:
        script = os.path.join(SRCDIR, "get_DNA_sequences.py")
    threads: 1
    resources:
        runtime = 60,
        mem_mb = MEM * 2
    log:
        "logs/get_DNA_sequences/{s}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "get_DNA_sequences: genes and GFF of the smORF candidates of {wildcards.s}"
    shell:
        """
        python {params.script} -p {input.p} -gd {input.gd} -gs {input.gs} -ad {input.ad} -asm {input.asm} \\
            -o {output.g} -a {output.a} > {log} 2>&1
        """
