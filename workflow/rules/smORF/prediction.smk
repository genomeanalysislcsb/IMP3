# smORF candidates of each sample (rules of the paper's smORF_prediction.smk), two ways (prediction: regions):
#   intergenic (default)  the CDS of the IMP3 annotation with <= max_aa aa, plus new genes called by pyrodigal-gv only
#                         on the intergenic regions of the IMP3 annotation (no second gene calling on the whole
#                         assembly); genes lifted back to contig coordinates
#   full                  the paper: pyrodigal-gv on the whole assembly
# Both: two pyrodigal-gv runs (default: genes >= 60 nt; small: genes >= 15 nt), smORF candidates of each (complete
# genes with an RBS motif, <= max_aa), deduplicated with MMseqs2 (100 % identity and coverage)

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


if INTERGENIC:
    rule intergenic_regions:
        input:
            asm = "Input/contigs/{s}.fasta.gz",
            gff = lambda wc: imp_file(wc.s, "annotation_gff")
        output:
            fasta = "Input/intergenic/{s}.fasta.gz",
            map = "Input/intergenic/{s}.map.tsv",
            stats = "Stats/{s}.intergenic.tsv"
        params:
            script = os.path.join(SRCDIR, "intergenic_regions.py"),
            occupied = PRED.get("occupied", "CDS tRNA tmRNA rRNA ncRNA CRISPR"),
            min_len = int(PRED.get("min_region", 30)),
            flank = int(PRED.get("flank", 0))
        threads: 1
        resources:
            runtime = 120,
            mem_mb = MEM * 4
        log:
            "logs/intergenic_regions/{s}.log"
        conda:
            os.path.join(ENVDIR, "IMP_smorf.yaml")
        message:
            "intergenic_regions: regions without IMP3 features on {wildcards.s}"
        shell:
            """
            python {params.script} regions --assembly {input.asm} --gff {input.gff} --occupied "{params.occupied}" \\
                --min-length {params.min_len} --flank {params.flank} --out-fasta {output.fasta} \\
                --out-map {output.map} --stats {output.stats} > {log} 2>&1
            """


    rule imp3_small_cds:
        input:
            asm = "Input/contigs/{s}.fasta.gz",
            gff = lambda wc: imp_file(wc.s, "annotation_gff"),
            faa = lambda wc: imp_file(wc.s, "bakta_faa")
        output:
            p = "pyrodigal_gv/imp3/{s}.smorf.faa",
            d = "pyrodigal_gv/imp3/{s}.fna",
            a = "pyrodigal_gv/imp3/{s}.gff"
        params:
            script = os.path.join(SRCDIR, "intergenic_regions.py"),
            max_aa = int(PRED.get("max_aa", 100)),
            use = "yes" if IMP3_CDS else "no"
        threads: 1
        resources:
            runtime = 120,
            mem_mb = MEM * 4
        log:
            "logs/imp3_small_cds/{s}.log"
        conda:
            os.path.join(ENVDIR, "IMP_smorf.yaml")
        message:
            "imp3_small_cds: CDS of the IMP3 annotation with <= {params.max_aa} aa on {wildcards.s}"
        shell:
            """
            if [[ {params.use} == yes ]]; then
                python {params.script} imp3 --sample {wildcards.s} --assembly {input.asm} --gff {input.gff} \\
                    --faa {input.faa} --max-aa {params.max_aa} --out-faa {output.p} --out-fna {output.d} \\
                    --out-gff {output.a} > {log} 2>&1
            else
                echo "imp3_cds: false" > {log}; : > {output.p}; : > {output.d}; echo "##gff-version  3" > {output.a}
            fi
            """


rule pyrodigal:
    input:
        seq = "Input/intergenic/{s}.fasta.gz" if INTERGENIC else "Input/contigs/{s}.fasta.gz",
        map = "Input/intergenic/{s}.map.tsv" if INTERGENIC else []
    output:
        a = "pyrodigal_gv/{run}/{s}.gff",
        d = "pyrodigal_gv/{run}/{s}.fna",
        p = "pyrodigal_gv/{run}/{s}.faa"
    wildcard_constraints:
        run = "default|small"
    params:
        mode = PRED.get("mode", "meta"),
        opts = lambda wc: PRED.get(wc.run, "--min-gene 60 --max-overlap 60" if wc.run == "default"
                                   else "--min-gene 15 --max-overlap 15"),
        lift = os.path.join(SRCDIR, "intergenic_regions.py") if INTERGENIC else "",
        raw = "pyrodigal_gv/{run}/{s}.fragments"
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
        if [[ -n "{params.lift}" ]]; then
            # intergenic regions: genes on the fragments, then lifted to contig coordinates (<contig>_<k>)
            pyrodigal-gv -i {input.seq} -p {params.mode} -d {params.raw}.fna -a {params.raw}.faa -f gff \\
                -o {params.raw}.gff -j {threads} -c -m --no-stop-codon {params.opts} > {log} 2>&1
            python {params.lift} lift --map {input.map} --faa-in {params.raw}.faa --faa-out {output.p} \\
                --fna-in {params.raw}.fna --fna-out {output.d} --gff-in {params.raw}.gff --gff-out {output.a} >> {log} 2>&1
            rm -f {params.raw}.fna {params.raw}.faa {params.raw}.gff
        else
            pyrodigal-gv -i {input.seq} -p {params.mode} -d {output.d} -a {output.p} -f gff -o {output.a} -j {threads} \\
                -c -m --no-stop-codon {params.opts} > {log} 2>&1
        fi
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
        small = "pyrodigal_gv/small/{s}.smorf.faa",
        imp3 = "pyrodigal_gv/imp3/{s}.smorf.faa" if INTERGENIC else []
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
        "deduplicate_smorfs: identical smORF candidates of {wildcards.s} (MMseqs2)"
    shell:
        """
        ins=""; for f in {input.default} {input.small} {input.imp3}; do [[ -s $f ]] && ins="$ins $f"; done
        mmseqs easy-cluster $ins {params.d} {params.tmp} -c 1.0 --min-seq-id 1.0 \\
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
        asm = "pyrodigal_gv/small/{s}.gff",
        gi = "pyrodigal_gv/imp3/{s}.fna" if INTERGENIC else [],
        ai = "pyrodigal_gv/imp3/{s}.gff" if INTERGENIC else []
    output:
        g = "potential_smORF/{s}.smorf.fna",
        a = "potential_smORF/{s}.smorf.gff"
    params:
        script = os.path.join(SRCDIR, "get_DNA_sequences.py"),
        extra = lambda wc, input: ("-gx %s -ax %s" % (input.gi, input.ai)) if INTERGENIC else ""
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
            -o {output.g} -a {output.a} {params.extra} > {log} 2>&1
        """
