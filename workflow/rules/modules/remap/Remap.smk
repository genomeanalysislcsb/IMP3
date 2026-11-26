# Remapping (steps: ... remap): the preprocessed reads of this run (e.g. an extracellular fraction) are mapped to the
# assembly of another IMP3 run of the same sample (e.g. the intracellular MG/MT co-assembly, remap: assembly_dir) and
# counted on its annotation and its MAGs, with the counting rules of the main workflow:
#   remapping/Assembly/<omic>.reads.sorted.bam                     reads on the other run's assembly (bwa-mem2)
#   remapping/Stats/<omic>/<omic>.{flagstat.txt,contig_coverage.tsv,strandedness.json,bin_abundance.tsv}
#   remapping/annotation/<omic>.{ID,Name}.<feature>_counts.tsv     genes (as Analysis/annotation)
#   remapping/annotation/<omic>.mantis.<db>.counts.tsv             terms of the merged annotation per database
#   remapping/annotation/per_contig_counts/...                     the same per contig
#   remapping/annotation/per_bin_counts/<omic>.mantis.<db>.per_bin.counts.tsv   Mantis terms per MAG
import glob as _glob

REMAP = config.get("remap", {}) or {}
REMAP_DIR = os.path.abspath(os.path.expandvars(str(REMAP.get("assembly_dir", ""))))
if not REMAP.get("assembly_dir") or not os.path.isdir(REMAP_DIR):
    raise Exception("steps: remap needs remap: assembly_dir, the output folder of the IMP3 run to map to (%s)" % REMAP_DIR)


def _remap_file(key, default, required=True):
    """a file of the other run, relative to remap: assembly_dir; .gz and plain both accepted"""
    rel = str(REMAP.get(key, default))
    for cand in [rel, rel[:-3] if rel.endswith(".gz") else rel + ".gz"]:
        hits = sorted(_glob.glob(os.path.join(REMAP_DIR, cand)))
        if hits:
            return hits[0]
    if required:
        raise Exception("remap: %s not found in %s (remap: %s)" % (rel, REMAP_DIR, key))
    return None


REMAP_ASSEMBLY = _remap_file("assembly", "Assembly/*.assembly.merged.fa.gz")
REMAP_GFF = _remap_file("gff", "Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz")
REMAP_BINS = _remap_file("bins", "Binning/MetaWrap/bins.contigs.gz", required=False)
REMAP_MANTIS_DBS = str(REMAP.get("mantis_dbs", "eggnog cazy enzyme_ec cog arcog go is_essential_gene kegg_brite kegg_ko "
                                 "kegg_module kegg_pathway kegg_reaction pfam tcdb tigrfam")).split()
REMAP_FEATURES = ["CRISPR", "CDS", "crispr-repeat", "crispr-spacer", "ncRNA", "oriC", "oriT", "regulatory_region",
                  "rRNA", "tmRNA", "tRNA"]
REMAP_INDEX = "remapping/Assembly/remap.assembly.fa"


def remap_mantis_gff(wildcards):
    """per-database GFF made from the other run's merged annotation (one line per term, Pfam names -> accessions)"""
    return "remapping/Assembly/terms/annotation.%s.gff" % wildcards.mdb


def remap_paired(wildcards):
    return (wildcards.type == "mg" and len(MG) >= 2) or (wildcards.type == "mt" and len(MT) >= 2)


def remap_reads(wildcards):
    reads = ["r1", "r2", "se"] if remap_paired(wildcards) else ["se"]
    return expand("Preprocessing/{type}.{read}.preprocessed.fq.gz", type=wildcards.type, read=reads)


def remap_stranding(wildcards):
    return str((REMAP.get("stranding") or {}).get(wildcards.type, config["featureCountsStranding"][wildcards.type])).lower()


def remap_fc_flags(wildcards):
    return "-p -O --fraction -M --fraction" if remap_paired(wildcards) else "-O --fraction -M --fraction"


REMAP_PFAM_HMM = (config.get("annotation_terms", {}) or {}).get("pfam_hmm", "") or \
    os.path.join(DBPATH, "mantis", "pfam", "Pfam-A.hmm")


rule remap_split_terms:
    input:
        gff = REMAP_GFF
    output:
        temp("remapping/Assembly/terms/annotation.{mdb}.gff")
    wildcard_constraints:
        mdb = "[A-Za-z0-9_-]+"
    params:
        script = os.path.join(SRCDIR, "split_annotation_terms.py"),
        pfam = ("--pfam-hmm %s" % REMAP_PFAM_HMM if os.path.exists(REMAP_PFAM_HMM) else "") + (" --pfam-metadata %s" % os.path.join(os.path.dirname(REMAP_PFAM_HMM), "metadata.tsv") if os.path.exists(os.path.join(os.path.dirname(REMAP_PFAM_HMM), "metadata.tsv")) else "")
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    log:
        "logs/remap/annotation_terms.{mdb}.log"
    message:
        "remap_split_terms: {wildcards.mdb} terms of the other run's merged annotation."
    shell:
        "python {params.script} --gff {input.gff} --db {wildcards.mdb} {params.pfam} --out {output} > {log} 2>&1"


rule remap_index:
    input:
        REMAP_ASSEMBLY
    output:
        expand(REMAP_INDEX + ".{ext}", ext=["0123", "amb", "ann", "bwt.2bit.64", "pac"])
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/remap/index.log"
    message:
        "remap_index: indexing the assembly of %s" % REMAP_DIR
    shell:
        """
        echo "Assembly: {input}" > {log}
        bwa-mem2 index -p {REMAP_INDEX} {input} >> {log} 2>&1
        """


rule remap_mapping:
    wildcard_constraints:
        type = "mg|mt"
    input:
        reads = remap_reads,
        index = rules.remap_index.output
    output:
        bam = "remapping/Assembly/{type}.reads.sorted.bam",
        bai = "remapping/Assembly/{type}.reads.sorted.bam.bai",
        flagstat = "remapping/Stats/{type}/{type}.flagstat.txt"
    params:
        paired = lambda wildcards: "yes" if remap_paired(wildcards) else "no",
        rg = lambda wildcards: "@RG\\tID:%s_%s_remap\\tSM:%s\\tLB:%s\\tPL:ILLUMINA" % (SAMPLE, wildcards.type, SAMPLE,
                                                                                 wildcards.type),
        se = lambda wildcards, input: input.reads[-1],
        mem = str(max(int(MAIN_MEMCORE_MB * 0.75), 1)) + "M"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/remap/mapping.{type}.log"
    message:
        "remap_mapping: {wildcards.type} reads on the assembly of %s" % REMAP_DIR
    shadow: "minimal"
    shell:
        """
        parts=""
        if [[ "{params.paired}" == "yes" ]]; then
            bwa-mem2 mem -v 1 -K {BWA_BATCH} -t {threads} -R '{params.rg}' {REMAP_INDEX} \\
                {input.reads[0]} {input.reads[1]} 2>> {log} | \\
                samtools sort -@ {threads} -m {params.mem} -o paired.bam - >> {log} 2>&1
            parts="paired.bam"
        fi
        se={params.se}
        if [[ -s $se ]]; then
            bwa-mem2 mem -v 1 -K {BWA_BATCH} -t {threads} -R '{params.rg}' {REMAP_INDEX} $se 2>> {log} | \\
                samtools sort -@ {threads} -m {params.mem} -o single.bam - >> {log} 2>&1
            parts="$parts single.bam"
        fi
        set -- $parts
        if [[ $# -gt 1 ]]; then
            samtools merge -@ {threads} -f {output.bam} $parts >> {log} 2>&1
        else
            mv $1 {output.bam}
        fi
        samtools index {output.bam} >> {log} 2>&1
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat} 2>> {log}
        """


rule remap_coverage:
    wildcard_constraints:
        type = "mg|mt"
    input:
        bam = "remapping/Assembly/{type}.reads.sorted.bam",
        bai = "remapping/Assembly/{type}.reads.sorted.bam.bai"
    output:
        "remapping/Stats/{type}/{type}.contig_coverage.tsv"
    resources:
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(4)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    log:
        "logs/remap/coverage.{type}.log"
    message:
        "remap_coverage: reads and depth per contig ({wildcards.type})."
    shell:
        "samtools coverage {input.bam} > {output} 2> {log}"


rule remap_strandedness:
    wildcard_constraints:
        type = "mg|mt"
    input:
        map = "remapping/Assembly/{type}.reads.sorted.bam",
        gff = REMAP_GFF
    output:
        txt = "remapping/Stats/{type}/{type}.strandedness.txt",
        json = "remapping/Stats/{type}/{type}.strandedness.json"
    params:
        stranding = remap_stranding,
        paired = lambda wildcards: "-p --countReadPairs" if remap_paired(wildcards) else "",
        min_ratio = 3
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(8)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/remap/strandedness.{type}.log"
    message:
        "remap_strandedness: strand of the {wildcards.type} library on the remapped reads."
    shadow: "minimal"
    shell:
        """
        if [[ "{params.stranding}" != "auto" ]]; then
            echo {params.stranding} > {output.txt}
            printf '{{"type": "%s", "stranding": %s, "detected": false}}\\n' {wildcards.type} {params.stranding} > {output.json}
            exit 0
        fi
        for s in 1 2; do
            featureCounts {params.paired} -t CDS -g ID -s $s -a {input.gff} -T {threads} -o fc.s$s.tsv {input.map} >> {log} 2>&1
        done
        s1=$(awk '$1 == "Assigned" {{print $2}}' fc.s1.tsv.summary)
        s2=$(awk '$1 == "Assigned" {{print $2}}' fc.s2.tsv.summary)
        choice=$(awk -v a="$s1" -v b="$s2" -v r={params.min_ratio} 'BEGIN {{ if (b >= r * a && b > 0) print 2; else if (a >= r * b && a > 0) print 1; else print 0 }}')
        echo $choice > {output.txt}
        printf '{{"type": "%s", "assigned_s1": %s, "assigned_s2": %s, "min_ratio": %s, "stranding": %s, "detected": true}}\\n' \\
            {wildcards.type} $s1 $s2 {params.min_ratio} $choice > {output.json}
        """


rule remap_fc_gene:
    input:
        map = "remapping/Assembly/{type}.reads.sorted.bam",
        strand = "remapping/Stats/{type}/{type}.strandedness.txt",
        gff = REMAP_GFF
    output:
        "remapping/annotation/{type}.{db}.{feature}_counts.tsv"
    wildcard_constraints:
        type = "mg|mt",
        db = "ID|Name"
    params:
        flags = remap_fc_flags,
        per_contig = ""
    resources:
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(12)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/remap/featurecount.{type}.{db}.{feature}.log"
    message:
        "remap_fc: {wildcards.type} reads, {wildcards.db} {wildcards.feature}."
    shadow: "minimal"
    shell:
        """
        featureCounts {params.flags} {params.per_contig} -t {wildcards.feature} -g {wildcards.db} \\
                      -s $(cat {input.strand}) -a {input.gff} -T {threads} -o {output} {input.map} > {log} 2>&1 \\
            || {{ echo "no {wildcards.feature} features counted" >> {log}; touch {output}; }}
        """


use rule remap_fc_gene as remap_fc_gene_per_contig with:
    output:
        "remapping/annotation/per_contig_counts/{type}.{db}.{feature}_per_contig.counts.tsv"
    params:
        flags = remap_fc_flags,
        per_contig = "-f"
    log:
        "logs/remap/featurecount.per_contig.{type}.{db}.{feature}.log"


rule remap_fc_mantis:
    wildcard_constraints:
        type = "mg|mt"
    input:
        map = "remapping/Assembly/{type}.reads.sorted.bam",
        strand = "remapping/Stats/{type}/{type}.strandedness.txt",
        gff = remap_mantis_gff
    output:
        "remapping/annotation/{type}.mantis.{mdb}.counts.tsv"
    params:
        flags = remap_fc_flags,
        per_contig = ""
    resources:
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(12)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/remap/featurecount.mantis.{type}.{mdb}.log"
    message:
        "remap_fc: {wildcards.type} reads, Mantis {wildcards.mdb}."
    shadow: "minimal"
    shell:
        """
        featureCounts {params.flags} {params.per_contig} -t CDS -g {wildcards.mdb} \\
                      -s $(cat {input.strand}) -a {input.gff} -T {threads} -o {output} {input.map} > {log} 2>&1 \\
            || {{ echo "no {wildcards.mdb} terms counted" >> {log}; touch {output}; }}
        """


use rule remap_fc_mantis as remap_fc_mantis_per_contig with:
    output:
        "remapping/annotation/per_contig_counts/{type}.mantis.{mdb}.per_contig.counts.tsv"
    params:
        flags = remap_fc_flags,
        per_contig = "-f"
    log:
        "logs/remap/featurecount.per_contig.mantis.{type}.{mdb}.log"


if REMAP_BINS:
    rule remap_bin_abundance:
        input:
            coverage = "remapping/Stats/{type}/{type}.contig_coverage.tsv",
            bins = REMAP_BINS
        output:
            "remapping/Stats/{type}/{type}.bin_abundance.tsv"
        params:
            script = os.path.join(SRCDIR, "remap_bins.py")
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            1
        log:
            "logs/remap/bin_abundance.{type}.log"
        message:
            "remap_bins: {wildcards.type} reads per MAG."
        shell:
            "python {params.script} abundance --coverage {input.coverage} --bins {input.bins} --out {output} > {log} 2>&1"

    rule remap_bin_functions:
        input:
            counts = "remapping/annotation/per_contig_counts/{type}.mantis.{mdb}.per_contig.counts.tsv",
            bins = REMAP_BINS
        output:
            "remapping/annotation/per_bin_counts/{type}.mantis.{mdb}.per_bin.counts.tsv"
        params:
            script = os.path.join(SRCDIR, "remap_bins.py")
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            1
        log:
            "logs/remap/bin_functions.{type}.{mdb}.log"
        message:
            "remap_bins: {wildcards.type} Mantis {wildcards.mdb} per MAG."
        shell:
            "python {params.script} functions --counts {input.counts} --bins {input.bins} --out {output} > {log} 2>&1"


rule REMAP:
    input:
        expand("remapping/Assembly/{type}.reads.sorted.bam.bai", type=TYPES),
        expand("remapping/Stats/{type}/{type}.{f}", type=TYPES,
               f=["flagstat.txt", "contig_coverage.tsv", "strandedness.json"]),
        expand("remapping/Stats/{type}/{type}.bin_abundance.tsv", type=TYPES) if REMAP_BINS else [],
        expand("remapping/annotation/{type}.{db}.{feature}_counts.tsv", type=TYPES, db=["ID", "Name"],
               feature=REMAP_FEATURES),
        expand("remapping/annotation/per_contig_counts/{type}.{db}.{feature}_per_contig.counts.tsv", type=TYPES,
               db=["ID", "Name"], feature=REMAP_FEATURES),
        expand("remapping/annotation/{type}.mantis.{mdb}.counts.tsv", type=TYPES, mdb=REMAP_MANTIS_DBS),
        expand("remapping/annotation/per_contig_counts/{type}.mantis.{mdb}.per_contig.counts.tsv", type=TYPES,
               mdb=REMAP_MANTIS_DBS),
        expand("remapping/annotation/per_bin_counts/{type}.mantis.{mdb}.per_bin.counts.tsv", type=TYPES,
               mdb=REMAP_MANTIS_DBS) if REMAP_BINS else []
    output:
        touch("status/remap.done")
