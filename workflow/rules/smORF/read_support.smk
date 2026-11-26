# metaG / metaT reads on the smORFs: featureCounts on the read alignments of the IMP3 run (reads on the same
# assembly; Assembly/<omic>.reads.sorted.bam), unstranded, multi-mapping and multi-overlapping reads as fractions
# (the settings of IMP3's gene counts). A smORF is supported by an omic with >= read_support_min_reads reads.

rule smorf_saf:
    input:
        "smORF/{s}.smorf.gff"
    output:
        "read_support/{s}.smorf.saf"
    params:
        script = os.path.join(SRCDIR, "smorf_read_support.py")
    threads: 1
    resources:
        runtime = 30,
        mem_mb = MEM
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        "python {params.script} saf --gff {input} --out {output} 2> /dev/null"


rule smorf_featurecounts:
    input:
        saf = "read_support/{s}.smorf.saf",
        bam = lambda wc: imp_file(wc.s, "bam", wc.omic)
    output:
        "read_support/{s}.{omic}.counts.tsv"
    wildcard_constraints:
        omic = "mg|mt"
    params:
        tmp = os.path.join(TMP, "fc_{s}_{omic}")
    threads: T_SMALL
    resources:
        runtime = 240,
        mem_mb = MEM * T_SMALL
    log:
        "logs/read_support/{s}.{omic}.log"
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "smorf_featurecounts: {wildcards.omic} reads of {wildcards.s} on its smORFs"
    shell:
        """
        mkdir -p {params.tmp}
        # paired-end reads: -p (as IMP3's gene counts); single-end otherwise
        n=$(samtools view -f 1 {input.bam} 2> /dev/null | head -n 1 | wc -l || true)
        if [[ $n -gt 0 ]]; then pe="-p"; else pe=""; fi
        featureCounts $pe -O --fraction -M --fraction -s 0 -F SAF -a {input.saf} -T {threads} --tmpDir {params.tmp} \\
            -o {output} {input.bam} > {log} 2>&1
        rm -rf {params.tmp}
        """


rule smorf_read_support:
    input:
        lambda wc: expand("read_support/{s}.{omic}.counts.tsv", s=wc.s, omic=OMICS[wc.s])
    output:
        tab = "read_support/{s}.smorf_reads.tsv",
        stats = "read_support/{s}.support_stats.tsv"
    params:
        script = os.path.join(SRCDIR, "smorf_read_support.py"),
        counts = lambda wc, input: " ".join("%s=%s" % (o, f) for o, f in zip(OMICS[wc.s], input)),
        min_reads = config.get("read_support_min_reads", 1)
    threads: 1
    resources:
        runtime = 30,
        mem_mb = MEM * 2
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        """
        python {params.script} summary --counts {params.counts} --min-reads {params.min_reads} \\
            --sample {wildcards.s} --out {output.tab} --stats {output.stats}
        """


rule smorf_read_support_samples:
    input:
        expand("read_support/{s}.support_stats.tsv", s=RS_SAMPLES)
    output:
        "read_support_samples.tsv"
    threads: 1
    resources:
        runtime = 10,
        mem_mb = MEM
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    shell:
        """
        python - {output} {input} <<'EOF'
import sys
out, files = sys.argv[1], sys.argv[2:]
cols, rows = ["sample", "omics", "smORFs", "MG only", "MT only", "MG+MT", "none"], []
for f in files:
    h, v = [l.rstrip("\\n").split("\\t") for l in open(f).read().splitlines()[:2]]
    d = dict(zip(h, v))
    rows.append([d.get(c, "0") for c in cols])
with open(out, "w") as fh:
    fh.write("\\t".join(cols) + "\\n")
    for r in rows:
        fh.write("\\t".join(r) + "\\n")
EOF
        """
