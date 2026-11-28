#!/usr/bin/env python3
"""smORFs of one IMP3 sample (imp3 smorf): one process, outputs in the sample's IMP3 folder.

Candidates
  intergenic (default)  the CDS of the IMP3 annotation with <= max_aa aa (complete, not touching a contig end) and
                        new genes called by pyrodigal-gv only in the regions without an IMP3 feature (types in
                        `occupied`), each region extended by `flank` nt; every region has a class:
                        between_features, contig_end (before the first / after the last feature of a contig),
                        unannotated_contig (contig without any feature)
  full                  pyrodigal-gv on the whole assembly, as in the smORF paper (Davin, Ortis Sunyer, Delgado et al.,
                        Nat Commun 2026); `paper: true` adds the paper's RBS-motif filter and its Rfam / Aragorn
                        filters (reproduction)
  gene calling          two pyrodigal-gv runs, metagenome mode, closed ends, masked N runs (as the paper):
                        default (min_gene 60, max_overlap 60) and small (15, 15); genes <= max_aa aa
Then
  exact duplicates      identical proteins merged (= MMseqs2 100 % identity and coverage); representative IMP3 >
                        default > small, then by position; `copies` = number of positions
  AntiFam               spurious ORFs (pyhmmer, gathering cut-offs) removed from the new candidates (Bakta checked
                        its own CDS)
  start status          start_confirmed: an in-frame stop codon upstream on the contig (the gene cannot be longer);
                        start_open: frame open to the contig end (may be the 3' part of a longer gene called from an
                        internal start codon; closed-end gene calling never reports partial genes). Core set =
                        start_confirmed
  RBS motif             reported, not filtered: motif / none / not_assessable (< min_upstream nt upstream in the
                        sequence the gene caller saw); IMP3 CDS: pyrodigal-gv on a window around the gene
  overlaps              IMP3 features overlapping a smORF (other than the smORF itself), e.g. riboswitch leaders
  metaP status          predicted -> validated with >= min_unique_peptides identified peptides found in the smORF and
                        in no protein > max_aa aa (IMP3 proteins, human, contaminants); I = L
  read support          featureCounts on the IMP3 BAMs (as IMP3's gene counts: -O --fraction -M --fraction, -s 0)
Outputs (in <sample>/Analysis/smORF/): <sample>.smorf.{faa,fna,gff}, <sample>.smorf_table.tsv,
<sample>.smorf_stats.tsv, <sample>.smorf.provenance.json, read counts; the IMP3 GFF with the smORFs as
<sample>/Analysis/annotation/annotation_CDS_RNA_hmms_smORF.gff.gz (the IMP3 GFF itself is not changed).
"""
import argparse
import bisect
import collections
import concurrent.futures as cf
import datetime
import glob
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time

VERSION = "1.0"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = {"TAA", "TAG", "TGA"}
CODE11 = dict(zip(
    [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
SET_RANK = {"imp3": 0, "default": 1, "small": 2}
T0 = time.time()


def log(msg):
    sys.stderr.write("[%7.1f s] %s\n" % (time.time() - T0, msg))
    sys.stderr.flush()


def zopen(path, mode="rt"):
    if "w" in mode:
        return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, mode) if gz else open(path, mode)


def first(*pats):
    for p in pats:
        for f in sorted(glob.glob(p)):
            if os.path.exists(f):
                return f
    return None


def read_fasta(path):
    name, seq = None, []
    with zopen(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(dna):
    """table 11, first codon as M (as pyrodigal / Bakta), stop codons removed at the end"""
    aa = "".join(CODE11.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))
    return ("M" + aa[1:]).rstrip("*") if aa else aa


def attrs(s):
    out = collections.OrderedDict()
    for x in s.strip().rstrip(";").split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def md5(path):
    if not path or not os.path.exists(path):
        return None
    if os.path.getsize(path) > 2 ** 31:          # large databases: size and time instead of a checksum
        return "size=%d;mtime=%d" % (os.path.getsize(path), os.path.getmtime(path))
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def tool_version(cmd):
    try:
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=30)
        return " ".join(r.stdout.split()[:12])
    except Exception:
        return None


# ------------------------------------------------------------------------------------------------- inputs
class Sample:
    def __init__(self, d, name):
        self.dir, self.name = d, name
        # the FASTA itself, not its indexes (.fa.0123, .fa.amb, .fa.fai, ... next to it)
        self.assembly = first(*[d + "/Assembly/%s.assembly.merged.%s" % (t, x) for t in ("mgmt", "mg", "mt")
                                for x in ("fa.gz", "fa")])
        self.gff = first(d + "/Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz",
                         d + "/Analysis/annotation/annotation_CDS_RNA_hmms.gff",
                         d + "/Analysis/annotation/bakta/annotation.filt.gff*")
        self.faa = first(d + "/Analysis/annotation/bakta/annotation.faa*")
        self.bams = {o: first(d + "/Assembly/%s.reads.sorted.bam" % o) for o in ("mg", "mt")}
        if not self.assembly or not self.gff:
            sys.exit("%s: assembly (%s) or IMP3 GFF (%s) not found" % (d, self.assembly, self.gff))


def read_gff(path):
    """header lines, [(contig, source, type, start, end, score, strand, phase, attr string)]"""
    head, rows = [], []
    with zopen(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                head.append(line.rstrip("\n"))
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) == 9:
                p[3], p[4] = int(p[3]), int(p[4])
                rows.append(p)
    return head, rows


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


# ------------------------------------------------------------------------------------------------- gene calling
def gene_finder(min_gene, max_overlap):
    import pyrodigal_gv
    kw = dict(meta=True, closed=True, mask=True)
    try:
        return pyrodigal_gv.ViralGeneFinder(min_gene=min_gene, max_overlap=max_overlap, **kw)
    except TypeError:
        sys.exit("pyrodigal-gv %s: no min_gene / max_overlap options (IMP_smorf.yaml pins 0.3.2)"
                 % getattr(pyrodigal_gv, "__version__", "?"))


def call_genes(finder, seqs, threads):
    """[(index, [genes])] for the sequences, in parallel (pyrodigal releases the GIL)"""
    with cf.ThreadPoolExecutor(max_workers=threads) as ex:
        return list(ex.map(finder.find_genes, seqs, chunksize=64))


def rbs_of(g):
    m = g.rbs_motif
    return (m if m and m != "None" else "None"), (g.rbs_spacer if g.rbs_spacer and g.rbs_spacer != "None" else "")


# ------------------------------------------------------------------------------------------------- main steps
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-dir", required=True, help="IMP3 output folder of the sample (imp/<sample>)")
    ap.add_argument("--sample", default="", help="sample name [folder name]")
    ap.add_argument("--config", default="", help="settings (YAML, config/config.smorf.yaml); options override")
    ap.add_argument("--db-path", default="", help="IMP3 databases (AntiFam / Rfam under <db>/smorf)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--outdir", default="", help="[<sample-dir>/Analysis/smORF]")
    a = ap.parse_args()
    # the environment's programs (featureCounts, samtools, cmscan, aragorn) next to its Python, also when the script
    # is started with that Python without activating the environment (imp3 smorf does so)
    os.environ["PATH"] = os.path.dirname(os.path.abspath(sys.executable)) + os.pathsep + os.environ.get("PATH", "")

    C = {"regions": "intergenic", "occupied": "CDS tRNA tmRNA rRNA ncRNA CRISPR", "min_region": 30, "flank": 0,
         "imp3_cds": True, "max_aa": 100, "default": [60, 60], "small": [15, 15], "antifam": True,
         "min_upstream": 20, "paper": False, "read_support": "auto", "read_support_min_reads": 1,
         "metap_peptides": "", "min_unique_peptides": 1, "metap_background": [], "databases": {},
         "pfam": "auto", "gmsc": "auto", "gmsc_mapper": "", "gmsc_tool": "diamond", "gmsc_min_identity": 0.9,
         "gmsc_min_coverage": 0.9, "gmsc_evalue": 1e-3}
    if a.config:
        C.update({k: v for k, v in read_config(a.config).items() if v is not None})
    S = Sample(os.path.abspath(a.sample_dir), a.sample or os.path.basename(os.path.abspath(a.sample_dir)))
    out = os.path.abspath(a.outdir or os.path.join(S.dir, "Analysis", "smORF"))
    os.makedirs(out, exist_ok=True)
    db = a.db_path or C.get("db_path", "")
    dbs = C.get("databases") or {}
    antifam_hmm = dbs.get("antifam_hmm") or os.path.join(db, "smorf", "AntiFam", "AntiFam.hmm")
    rfam_cm = dbs.get("rfam_cm") or os.path.join(db, "smorf", "Rfam", "Rfam.cm")
    full = C["regions"] == "full"
    paper = bool(C.get("paper")) and full
    max_aa, mu = int(C["max_aa"]), int(C["min_upstream"])
    stats = collections.OrderedDict(sample=S.name, regions=C["regions"], flank=int(C["flank"]))
    log("%s: %s mode, %d threads -> %s" % (S.name, C["regions"], a.threads, out))

    # assembly and IMP3 annotation
    seq = {n: s.upper() for n, s in read_fasta(S.assembly)}
    order = list(seq)
    stats["contigs"], stats["assembly_bp"] = len(seq), sum(len(s) for s in seq.values())
    gff_head, gff_rows = read_gff(S.gff)
    occupied = set(str(C["occupied"]).split())
    feats = collections.defaultdict(list)          # contig -> sorted [(start, end, type, label)]
    occ = collections.defaultdict(list)
    for r in gff_rows:
        c, t, s0, e0, at = r[0], r[2], r[3], r[4], r[8]
        if c not in seq or t == "region":
            continue
        d = attrs(at)
        feats[c].append((s0, e0, t, (d.get("product") or d.get("Name") or d.get("ID") or "").replace(" ", "_")))
        if t in occupied:
            occ[c].append((s0, e0))
    for c in feats:
        feats[c].sort()
    log("assembly %d contigs, %d bp; IMP3 features %d" % (len(seq), stats["assembly_bp"], len(gff_rows)))

    # regions (contig, frag_start, frag_end, core_start, core_end, class)
    if full:
        frags = [(c, 1, len(s), 1, len(s), "whole_contig") for c, s in seq.items()]
    else:
        frags, fl = [], int(C["flank"])
        for c, s in seq.items():
            L, pos, gaps = len(s), 1, []
            mo = merge(occ.get(c, []))
            for s0, e0 in mo + [[L + 1, L + 1]]:
                if s0 - pos >= int(C["min_region"]):
                    gaps.append((pos, s0 - 1))
                pos = max(pos, e0 + 1)
            for s0, e0 in gaps:
                cls = ("unannotated_contig" if not mo else
                       "contig_end" if (s0 == 1 or e0 == L) else "between_features")
                frags.append((c, max(1, s0 - fl), min(L, e0 + fl), s0, e0, cls))
        rc = collections.Counter(f[5] for f in frags)
        rb = collections.Counter()
        for f in frags:
            rb[f[5]] += f[4] - f[3] + 1
        stats["regions"] = len(frags)
        stats["regions_bp"] = sum(rb.values())
        for k in ("between_features", "contig_end", "unannotated_contig"):
            stats["regions_%s" % k], stats["regions_%s_bp" % k] = rc[k], rb[k]
        log("regions: %d (%s)" % (len(frags), ", ".join("%s %d / %.1f Mb" % (k, rc[k], rb[k] / 1e6) for k in rc)))

    # candidates: dict per smORF
    cands = []
    for setname, (mg, mo_) in (("default", C["default"]), ("small", C["small"])):
        fdr = gene_finder(int(mg), int(mo_))
        fseqs = [seq[f[0]][f[1] - 1:f[2]] for f in frags]
        res = call_genes(fdr, fseqs, a.threads)
        n = 0
        for f, genes in zip(frags, res):
            for g in genes:
                prot = g.translate().rstrip("*")
                if len(prot) > max_aa or g.partial_begin or g.partial_end:
                    continue
                s0, e0 = f[1] + g.begin - 1, f[1] + g.end - 1
                st = "+" if g.strand == 1 else "-"
                up_seen = (g.begin - 1) if st == "+" else (f[2] - f[1] + 1 - g.end)
                motif, spacer = rbs_of(g)
                cands.append(dict(set=setname, contig=f[0], start=s0, end=e0, strand=st, prot=prot,
                                  region=f[5], rbs_motif=motif, rbs_spacer=spacer, up_seen=up_seen,
                                  locus_tag="", gff_id=""))
                n += 1
        stats["candidates_%s" % setname] = n
        log("pyrodigal-gv %s (min_gene %s, max_overlap %s): %d candidates <= %d aa" % (setname, mg, mo_, n, max_aa))
    if not full and C.get("imp3_cds", True):
        n = 0
        small_finder = gene_finder(15, 15)
        wins, idx = [], []
        for r in gff_rows:
            c, t, s0, e0, st = r[0], r[2], r[3], r[4], r[6]
            if t != "CDS" or c not in seq:
                continue
            L = len(seq[c])
            at = r[8].lower()
            if "truncated" in at or re.search(r"partial=(1|01|10|11)", at) or s0 <= 3 or e0 >= L - 2:
                continue
            dna = seq[c][s0 - 1:e0]
            dna = revcomp(dna) if st == "-" else dna
            prot = translate(dna)
            if not prot or len(prot) > max_aa or (e0 - s0 + 1) % 3:
                continue
            d = attrs(r[8])
            cands.append(dict(set="imp3", contig=c, start=s0, end=e0, strand=st, prot=prot, region="imp3_cds",
                              rbs_motif="", rbs_spacer="", up_seen=0, locus_tag=d.get("locus_tag", ""),
                              gff_id=d.get("ID", "")))
            # window for the RBS motif: the gene plus 60 nt upstream (pyrodigal-gv, small settings)
            ws, we = (max(1, s0 - 60), e0) if st == "+" else (s0, min(L, e0 + 60))
            wins.append(seq[c][ws - 1:we])
            idx.append((len(cands) - 1, ws, we))
            n += 1
        stats["candidates_imp3"] = n
        res = call_genes(small_finder, wins, a.threads)
        matched = 0
        for (i, ws, we), genes in zip(idx, res):
            x = cands[i]
            for g in genes:
                gs, ge = ws + g.begin - 1, ws + g.end - 1
                if (gs, ge, "+" if g.strand == 1 else "-") == (x["start"], x["end"], x["strand"]):
                    x["rbs_motif"], x["rbs_spacer"] = rbs_of(g)
                    x["up_seen"] = (g.begin - 1) if x["strand"] == "+" else (we - ws + 1 - g.end)
                    matched += 1
                    break
            else:
                x["rbs_motif"] = "unmatched"
        stats["imp3_cds_rbs_matched"] = matched
        log("IMP3 CDS <= %d aa: %d (RBS from the same gene call: %d)" % (max_aa, n, matched))
    stats["candidates"] = len(cands)

    # exact duplicates
    by_prot = collections.defaultdict(list)
    for x in cands:
        by_prot[x["prot"]].append(x)
    smorfs = []
    for prot, xs in by_prot.items():
        xs.sort(key=lambda x: (SET_RANK[x["set"]], x["contig"], x["start"]))
        rep = xs[0]
        rep["copies"] = len({(x["contig"], x["start"], x["end"], x["strand"]) for x in xs})
        rep["found_by"] = ",".join(sorted({x["set"] for x in xs}, key=SET_RANK.get))
        smorfs.append(rep)
    stats["deduplicated"] = len(smorfs)
    log("exact duplicates merged: %d smORF candidates" % len(smorfs))

    # filters
    removed = collections.Counter()
    if C.get("antifam", True):
        import pyhmmer
        with pyhmmer.plan7.HMMFile(antifam_hmm) as fh:
            hmms = list(fh)
        alph = hmms[0].alphabet
        test = [i for i, x in enumerate(smorfs) if x["set"] != "imp3" or full]
        if not test:
            test = []
        block = pyhmmer.easel.DigitalSequenceBlock(alph, [
            pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=smorfs[i]["prot"]).digitize(alph)
            for i in test])
        bad = set()
        for hits in (pyhmmer.hmmsearch(hmms, block, cpus=a.threads, bit_cutoffs="gathering") if test else []):
            for h in hits:
                if getattr(h, "reported", getattr(h, "is_reported", True)):
                    bad.add(int(h.name.decode() if isinstance(h.name, bytes) else h.name))
        for i in bad:
            smorfs[i]["drop"] = "antifam"
        removed["antifam"] = len(bad)
        log("AntiFam (%d HMMs, gathering): %d spurious of %d tested" % (len(hmms), len(bad), len(test)))
    if paper:
        for x in smorfs:
            if x["rbs_motif"] == "None" and "drop" not in x:
                x["drop"] = "rbs"
                removed["rbs"] += 1
        tmp = tempfile.mkdtemp(prefix="smorf_", dir=out)
        fna = os.path.join(tmp, "c.fna")
        with open(fna, "w") as fh:
            for i, x in enumerate(smorfs):
                if "drop" not in x:
                    fh.write(">s%d\n%s\n" % (i, dna_of(seq, x)))
        tbl = os.path.join(tmp, "rfam.tbl")
        subprocess.run(["cmscan", "--rfam", "--fmt", "2", "--nohmmonly", "--cut_tc", "-g", "--cpu", str(a.threads),
                        "--noali", "--tblout", tbl, "-o", os.devnull, rfam_cm, fna], check=True)
        hits = {l.split()[3] for l in open(tbl) if not l.startswith("#") and len(l.split()) > 3}
        for h in hits:
            smorfs[int(h[1:])]["drop"] = "rfam"
        removed["rfam"] = len(hits)
        ara = set()
        for code in ("-gcstd", "-gcbact"):
            o = os.path.join(tmp, "aragorn%s.txt" % code)
            subprocess.run(["aragorn", code, "-d", "-w", "-fasta", "-o", o, fna], check=True)
            lines = open(o).read().splitlines()
            for j, l in enumerate(lines[:-1]):
                if l.startswith(">") and "0 genes found" not in lines[j + 1]:
                    ara.add(l[1:].split()[0])
        for h in ara:
            if "drop" not in smorfs[int(h[1:])]:
                smorfs[int(h[1:])]["drop"] = "aragorn"
                removed["aragorn"] += 1
        shutil.rmtree(tmp, ignore_errors=True)
        log("paper filters: RBS %d, Rfam %d, Aragorn %d" % (removed["rbs"], removed["rfam"], removed["aragorn"]))
    for k in ("antifam", "rbs", "rfam", "aragorn"):
        stats["removed_%s" % k] = removed[k]
    smorfs = [x for x in smorfs if "drop" not in x]
    smorfs.sort(key=lambda x: (x["contig"], x["start"], x["strand"]))

    # start status, RBS status, region class in full mode, overlaps
    cds_iv = {c: sorted((f[0], f[1]) for f in fs if f[2] == "CDS") for c, fs in feats.items()}
    for k, x in enumerate(smorfs, 1):
        x["id"] = "%s_smORF_%d" % (S.name, k)
        c, s0, e0, st = x["contig"], x["start"], x["end"], x["strand"]
        x["start_status"] = "confirmed" if upstream_stop(seq[c], s0, e0, st) else "open"
        x["rbs_status"] = ("not_assessable" if x["rbs_motif"] == "unmatched" or x["up_seen"] < mu else
                           "none" if x["rbs_motif"] == "None" else "motif")
        if full:
            x["region"] = where(c, s0, e0, cds_iv.get(c, []), len(seq[c]), bool(occ.get(c)))
        ov = []
        fs = feats.get(c, [])
        j = bisect.bisect_left(fs, (s0 - 100000,))
        for f in fs[max(0, j - 1):]:
            if f[0] > e0:
                break
            if f[1] >= s0 and not (f[0] == s0 and f[1] == e0 and f[2] == "CDS"):
                ov.append("%s:%s" % (f[2], f[3] or "-"))
        x["overlaps"] = ",".join(sorted(set(ov)))
    stats["smorfs"] = len(smorfs)
    for k, v in collections.Counter(x["set"] for x in smorfs).items():
        stats["smorfs_%s" % k] = v
    for k, v in collections.Counter(x["start_status"] for x in smorfs).items():
        stats["start_%s" % k] = v
    for k, v in collections.Counter(x["rbs_status"] for x in smorfs).items():
        stats["rbs_%s" % k] = v
    for k, v in collections.Counter(x["region"] for x in smorfs if x["set"] != "imp3").items():
        stats["new_in_%s" % k] = v
    L = sorted(len(x["prot"]) for x in smorfs)
    if L:
        stats.update(length_min=L[0], length_median=L[len(L) // 2], length_max=L[-1],
                     smorfs_le50aa=sum(v <= 50 for v in L), smorfs_lt20aa=sum(v < 20 for v in L))
    log("smORFs: %d (start confirmed %d, open %d)" % (len(smorfs), stats.get("start_confirmed", 0),
                                                      stats.get("start_open", 0)))

    # evidence from known proteins: Pfam domains (new smORFs) and GMSC (all smORFs); columns, never filters
    pfam_hmm = pfam_db(C, db)
    pfam(smorfs, pfam_hmm, C, a.threads, stats)
    gmsc_used = gmsc(smorfs, C, db, out, a.threads, stats)

    # metaP status
    peps_file = C.get("metap_peptides") or ""
    if peps_file and "{sample}" in peps_file:
        peps_file = peps_file.replace("{sample}", S.name)
    metap(S, smorfs, peps_file, C, max_aa, stats)

    # read support
    rs = str(C.get("read_support", "auto")).lower()
    bams = {o: b for o, b in S.bams.items() if b}
    if rs in ("true", "yes") and not bams:
        sys.exit("read_support: true, but no BAM in %s/Assembly" % S.dir)
    if rs not in ("false", "no") and bams:
        read_support(S, smorfs, bams, out, a.threads, float(C.get("read_support_min_reads", 1)), stats)

    # outputs
    # the GFF with smORFs goes next to the IMP3 GFF only for the standard output folder (not for test runs)
    gff_dst = (os.path.join(S.dir, "Analysis", "annotation", "annotation_CDS_RNA_hmms_smORF.gff.gz") if not a.outdir
               else os.path.join(out, "%s.annotation_CDS_RNA_hmms_smORF.gff.gz" % S.name))
    write_outputs(S, smorfs, seq, gff_head, gff_rows, out, stats, gff_dst)
    prov = collections.OrderedDict(
        tool="imp3_smorf.py", version=VERSION, date=datetime.datetime.now().isoformat(timespec="seconds"),
        host=platform.node(), python=platform.python_version(), threads=a.threads, settings=C,
        inputs={k: {"path": p, "bytes": os.path.getsize(p), "mtime": os.path.getmtime(p)}
                for k, p in (("assembly", S.assembly), ("gff", S.gff), ("bakta_faa", S.faa),
                             ("bam_mg", S.bams["mg"]), ("bam_mt", S.bams["mt"])) if p},
        databases={"antifam_hmm": {"path": antifam_hmm, "md5": md5(antifam_hmm)} if C.get("antifam", True) else None,
                   "rfam_cm": {"path": rfam_cm, "md5": md5(rfam_cm)} if paper else None,
                   "pfam_hmm": {"path": pfam_hmm, "md5": md5(pfam_hmm)} if pfam_hmm else None,
                   "gmsc": gmsc_used},
        software={"pyrodigal_gv": mod_version("pyrodigal_gv"), "pyrodigal": mod_version("pyrodigal"),
                  "pyhmmer": mod_version("pyhmmer"), "featureCounts": tool_version(["featureCounts", "-v"]),
                  "cmscan": tool_version(["cmscan", "-h"]) if paper else None,
                  "aragorn": tool_version(["aragorn", "-h"]) if paper else None},
        script_md5=md5(os.path.abspath(__file__)), runtime_s=round(time.time() - T0, 1), stats=stats)
    with open(os.path.join(out, "%s.smorf.provenance.json" % S.name), "w") as fh:
        json.dump(prov, fh, indent=1, default=str)
    log("done")


def _value(v):
    v = v.strip()
    if v.startswith("[") and v.endswith("]"):
        return [_value(x) for x in v[1:-1].split(",") if x.strip()]
    if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
        return v[1:-1]
    if v.lower() in ("true", "false"):
        return v.lower() == "true"
    if v in ("", "null", "~"):
        return None
    for t in (int, float):
        try:
            return t(v)
        except ValueError:
            pass
    return v


def read_config(path):
    """the settings file (YAML); PyYAML if installed, else the subset the settings use (key: value, [lists],
    one level of nesting, # comments)"""
    try:
        import yaml
        return yaml.safe_load(open(path)) or {}
    except ImportError:
        pass
    out, cur = {}, None
    for line in open(path):
        txt = re.sub(r"\s+#.*$", "", line.rstrip("\n")) if not line.lstrip().startswith("#") else ""
        if not txt.strip():
            continue
        m = re.match(r"^(\s*)([\w.-]+):\s*(.*)$", txt)
        if not m:
            continue
        ind, k, v = m.groups()
        if ind and cur is not None:
            out[cur][k] = _value(v)
        elif v.strip() == "":
            cur, out[k] = k, {}
        else:
            cur, out[k] = None, _value(v)
    return out


def mod_version(m):
    try:
        return __import__(m).__version__
    except Exception:
        return None


def dna_of(seq, x):
    d = seq[x["contig"]][x["start"] - 1:x["end"]]
    return revcomp(d) if x["strand"] == "-" else d


def upstream_stop(s, s0, e0, strand):
    if strand == "+":
        i = s0 - 4
        while i >= 0:
            if s[i:i + 3] in STOPS:
                return True
            i -= 3
        return False
    i = e0
    while i + 3 <= len(s):
        if revcomp(s[i:i + 3]) in STOPS:
            return True
        i += 3
    return False


def where(c, s0, e0, cds, L, annotated):
    j = bisect.bisect_right(cds, (e0, 10 ** 12))
    for a, b in cds[max(0, j - 50):j]:
        if a <= s0 and e0 <= b:
            return "inside_IMP3_CDS"
    for a, b in cds[max(0, j - 50):j]:
        if s0 <= b and a <= e0:
            return "overlaps_IMP3_CDS"
    if not annotated:
        return "unannotated_contig"
    before = any(b < s0 for a, b in cds[:j])
    after = j < len(cds)
    return "between_features" if before and after else "contig_end"


def pfam_db(C, db):
    """Pfam-A HMMs: databases.pfam_hmm, else <db>/smorf/Pfam/Pfam-A.hmm, else the Pfam-A of the Mantis databases"""
    if str(C.get("pfam", "auto")).lower() in ("false", "no"):
        return None
    p = (C.get("databases") or {}).get("pfam_hmm") or ""
    cands = [p] if p else [os.path.join(db, "smorf", "Pfam", "Pfam-A.hmm")] + sorted(
        glob.glob(os.path.join(db, "mantis", "*pfam*", "Pfam-A.hmm")) +
        glob.glob(os.path.join(os.path.dirname(db.rstrip("/")), "mantis", "*", "pfam", "Pfam-A.hmm")))
    for c in cands:
        if c and os.path.exists(c):
            return c
    if str(C.get("pfam")).lower() in ("true", "yes"):
        sys.exit("pfam: true, but no Pfam-A.hmm (set databases: pfam_hmm)")
    log("Pfam: no Pfam-A.hmm found (databases: pfam_hmm), skipped")
    return None


def pfam(smorfs, hmm, C, threads, stats):
    """Pfam-A domains (gathering cut-offs) of the new smORFs; IMP3 CDS keep their Mantis annotation in the GFF"""
    for x in smorfs:
        x["pfam"] = ""
    if not hmm:
        return
    import pyhmmer
    test = [i for i, x in enumerate(smorfs) if x["set"] != "imp3"]
    if not test:
        return
    with pyhmmer.plan7.HMMFile(hmm) as fh:
        hmms = list(fh)
    alph = hmms[0].alphabet
    block = pyhmmer.easel.DigitalSequenceBlock(alph, [
        pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=smorfs[i]["prot"]).digitize(alph) for i in test])
    hits = collections.defaultdict(set)
    for th in pyhmmer.hmmsearch(hmms, block, cpus=threads, bit_cutoffs="gathering"):
        q = th.query
        acc = getattr(q, "accession", None) or b""
        nm = getattr(q, "name", None) or b""
        acc, nm = (acc.decode() if isinstance(acc, bytes) else acc), (nm.decode() if isinstance(nm, bytes) else nm)
        for h in th:
            if getattr(h, "included", getattr(h, "is_included", True)):
                hits[int(h.name.decode() if isinstance(h.name, bytes) else h.name)].add("%s:%s" % (acc.split(".")[0], nm))
    for i, v in hits.items():
        smorfs[i]["pfam"] = ",".join(sorted(v))
    stats["pfam_hits"] = len(hits)
    log("Pfam (%d HMMs, gathering): %d of %d new smORFs with a domain" % (len(hmms), len(hits), len(test)))


def gmsc_mapper_cmd(C, db):
    """the gmsc-mapper program: gmsc_mapper, else <db>/gmsc/gmsc-mapper.path (written by the installer), else PATH"""
    p = C.get("gmsc_mapper") or ""
    f = os.path.join(db, "gmsc", "gmsc-mapper.path")
    if not p and os.path.exists(f):
        p = open(f).read().strip()
    return p if p and os.path.exists(p) else shutil.which("gmsc-mapper")


def gmsc(smorfs, C, db, out, threads, stats):
    """GMSC (global microbial smORF catalogue, Duan et al. 2024) with GMSC-mapper on the smORF proteins: best hit,
    'family' with >= gmsc_min_identity identity (GMSC's family definition, 90 %), else 'homolog'; habitat,
    taxonomy, quality label and CDD of the hits as GMSC-mapper reports them"""
    for x in smorfs:
        x.update(gmsc_class="", gmsc_family="", gmsc_identity="", gmsc_coverage="", gmsc_quality="", gmsc_habitat="",
                 gmsc_taxonomy="", gmsc_cdd="")
    mode = str(C.get("gmsc", "auto")).lower()
    if mode in ("false", "no"):
        return None
    dbdir = (C.get("databases") or {}).get("gmsc_db") or os.path.join(db, "gmsc")
    prog = gmsc_mapper_cmd(C, db)
    if not prog or not os.path.exists(os.path.join(dbdir, "GMSC10.90AA.faa.gz")):
        if mode in ("true", "yes"):
            sys.exit("gmsc: true, but gmsc-mapper (%s) or %s/GMSC10.90AA.faa.gz is missing" % (prog, dbdir))
        log("GMSC: gmsc-mapper or its database (%s) not found, skipped" % dbdir)
        return None
    tmp = tempfile.mkdtemp(prefix="gmsc_", dir=out)
    faa = os.path.join(tmp, "smorfs.faa")
    with open(faa, "w") as fh:
        for x in smorfs:
            fh.write(">%s\n%s\n" % (x["id"], x["prot"]))
    res = os.path.join(tmp, "out")
    cmd = [prog, "--aa-genes", faa, "-o", res, "--dbdir", dbdir, "--tool", str(C.get("gmsc_tool", "diamond")),
           "-t", str(threads), "--id", "0", "--cov", str(C.get("gmsc_min_coverage", 0.9)),
           "-e", str(C.get("gmsc_evalue", 1e-3)), "--quiet"]
    # gmsc-mapper calls diamond / mmseqs: its environment's bin first on PATH
    env = dict(os.environ, PATH=os.path.dirname(os.path.realpath(prog)) + os.pathsep + os.environ.get("PATH", ""))
    with open(os.path.join(out, "gmsc-mapper.log"), "w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env)
    if r.returncode:
        shutil.rmtree(tmp, ignore_errors=True)
        sys.exit("gmsc-mapper failed (exit %d), see %s/gmsc-mapper.log" % (r.returncode, out))

    def table(name):
        f = os.path.join(res, name)
        rows = []
        if os.path.exists(f):
            hdr = None
            for l in open(f):
                if l.startswith("#") or not l.strip():
                    continue
                p = l.rstrip("\n").split("\t")
                if hdr is None and p[0] == "qseqid":
                    hdr = p
                    continue
                rows.append(dict(zip(hdr or [], p)) if hdr else p)
        return rows

    best = {}
    for r in table("alignment.out.smorfs.tsv"):
        if not isinstance(r, dict) or "qseqid" not in r:
            continue
        b = best.get(r["qseqid"])
        if b is None or float(r.get("bitscore", 0)) > float(b.get("bitscore", 0)):
            best[r["qseqid"]] = r
    ann = {}
    for name, col in (("habitat.out.smorfs.tsv", "habitat"), ("taxonomy.out.smorfs.tsv", "taxonomy"),
                      ("quality.out.smorfs.tsv", "quality"), ("domain.out.smorfs.tsv", "cdd")):
        for r in table(name):
            if isinstance(r, dict):
                ann.setdefault(r.get("qseqid"), {})[col] = r.get(col, "")
            elif len(r) >= 2:
                ann.setdefault(r[0], {})[col] = r[1]
    mi = 100 * float(C.get("gmsc_min_identity", 0.9))
    for x in smorfs:
        b = best.get(x["id"])
        if not b:
            continue
        pid = float(b.get("pident", 0))
        x.update(gmsc_class="family" if pid >= mi else "homolog", gmsc_family=b.get("sseqid", ""),
                 gmsc_identity="%.1f" % pid, gmsc_coverage=b.get("qcovhsp", ""))
        for col in ("habitat", "taxonomy", "quality", "cdd"):
            x["gmsc_" + col] = (ann.get(x["id"]) or {}).get(col, "")
    c = collections.Counter(x["gmsc_class"] for x in smorfs)
    stats["gmsc_family"], stats["gmsc_homolog"] = c["family"], c["homolog"]
    stats["gmsc_high_quality"] = sum("high" in x["gmsc_quality"].lower() for x in smorfs if x["gmsc_class"])
    shutil.rmtree(tmp, ignore_errors=True)
    log("GMSC: %d smORFs in a GMSC family (>= %.0f %% identity), %d homologs, %d with a high-quality hit"
        % (c["family"], mi, c["homolog"], stats["gmsc_high_quality"]))
    return {"dbdir": dbdir, "faa": md5(os.path.join(dbdir, "GMSC10.90AA.faa.gz")),
            "gmsc_mapper": tool_version([prog, "--version"]), "tool": C.get("gmsc_tool", "diamond")}


def metap(S, smorfs, peps_file, C, max_aa, stats):
    for x in smorfs:
        x["status"], x["unique_peptides"], x["peptides"] = "predicted", 0, ""
    if not peps_file:
        stats["metap_measured"] = "no"
        return
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from smorf_metap import read_peptides, norm, K, read_fasta as rf
    peps = read_peptides(peps_file, S.name)
    by_pre = collections.defaultdict(list)
    for p in peps:
        by_pre[norm(p)[:K]].append(norm(p))

    def scan(s):
        s, hits = norm(s), set()
        for i in range(len(s) - K + 1):
            for p in by_pre.get(s[i:i + K], ()):
                if s.startswith(p, i):
                    hits.add(p)
        return hits

    hits = [scan(x["prot"]) for x in smorfs]
    found = set().union(*hits) if hits else set()
    shared = set()
    if found:
        by_pre = collections.defaultdict(list)
        for p in found:
            by_pre[p[:K]].append(p)
        for f in [S.faa] + list(C.get("metap_background") or []):
            if f:
                for _, s in rf(f):
                    if len(s.replace("*", "")) > max_aa:
                        shared |= scan(s)
    nv = 0
    for x, h in zip(smorfs, hits):
        u = sorted(h - shared)
        x["unique_peptides"], x["peptides"] = len(u), ",".join(u)
        if len(u) >= int(C.get("min_unique_peptides", 1)) and u:
            x["status"] = "validated"
            nv += 1
    stats.update(metap_measured="yes", peptides=len(peps), peptides_in_smorfs=len(found),
                 peptides_also_in_background=len(shared), validated=nv)
    log("metaP: %d peptides, %d smORFs validated" % (len(peps), nv))


def read_support(S, smorfs, bams, out, threads, min_reads, stats):
    tmp = tempfile.mkdtemp(prefix="fc_", dir=out)
    saf = os.path.join(tmp, "smorf.saf")
    with open(saf, "w") as fh:
        fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for x in smorfs:
            fh.write("%s\t%s\t%d\t%d\t%s\n" % (x["id"], x["contig"], x["start"], x["end"], x["strand"]))
    for o, bam in bams.items():
        pe = subprocess.run("samtools view -f 1 %s 2> /dev/null | head -n 1 | wc -l" % bam, shell=True,
                            stdout=subprocess.PIPE, universal_newlines=True).stdout.strip() not in ("", "0")
        res = os.path.join(out, "%s.smorf.%s.counts.tsv" % (S.name, o))
        cmd = ["featureCounts"] + (["-p"] if pe else []) + ["-O", "--fraction", "-M", "--fraction", "-s", "0", "-F",
               "SAF", "-a", saf, "-T", str(threads), "--tmpDir", tmp, "-o", res, bam]
        with open(os.path.join(tmp, "fc.%s.log" % o), "w") as lf:
            subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, check=True)
        cnt = {}
        for l in open(res):
            if l.startswith("#") or l.startswith("Geneid"):
                continue
            p = l.rstrip("\n").split("\t")
            cnt[p[0]] = float(p[-1])
        for x in smorfs:
            x["%s_reads" % o] = cnt.get(x["id"], 0.0)
        log("featureCounts %s (%s): %d smORFs with >= %g reads" % (o, "paired" if pe else "single",
                                                                   sum(v >= min_reads for v in cnt.values()), min_reads))
    for x in smorfs:
        sup = [o.upper() for o in ("mg", "mt") if x.get("%s_reads" % o, 0) >= min_reads and o in bams]
        x["support"] = ("MG+MT" if len(sup) == 2 else (sup[0] + " only") if sup else "none")
    for k, v in collections.Counter(x["support"] for x in smorfs).items():
        stats["support_%s" % k.replace(" ", "_").replace("+", "_")] = v
    shutil.rmtree(tmp, ignore_errors=True)


def write_outputs(S, smorfs, seq, gff_head, gff_rows, out, stats, dst):
    n, p = S.name, os.path.join(out, S.name)
    cols = ["smORF", "sample", "set", "found_by", "copies", "contig", "start", "end", "strand", "length_aa",
            "region", "start_status", "rbs_status", "rbs_motif", "rbs_spacer", "imp3_locus_tag", "status",
            "unique_peptides", "overlaps", "pfam", "gmsc_class", "gmsc_family", "gmsc_identity", "gmsc_coverage",
            "gmsc_quality", "gmsc_habitat", "gmsc_taxonomy", "gmsc_cdd", "mg_reads", "mt_reads", "support"]

    def gattr(x):
        a = collections.OrderedDict(
            smorf_id=x["id"], smorf_set=x["set"], smorf_length_aa=len(x["prot"]), smorf_region=x["region"],
            smorf_start=x["start_status"], smorf_rbs=x["rbs_status"], smorf_rbs_motif=x["rbs_motif"] or "-",
            smorf_status=x["status"], smorf_unique_peptides=x["unique_peptides"])
        if x["overlaps"]:
            a["smorf_overlaps"] = x["overlaps"]
        if x.get("pfam"):
            a["smorf_pfam"] = x["pfam"]
        if x.get("gmsc_class"):
            a["smorf_gmsc"] = x["gmsc_class"]
            a["smorf_gmsc_family"] = x["gmsc_family"]
            if x.get("gmsc_quality"):
                a["smorf_gmsc_quality"] = x["gmsc_quality"].replace(" ", "_").replace(";", ",")
        return ";".join("%s=%s" % kv for kv in a.items())

    with open(p + ".smorf.faa", "w") as fa, open(p + ".smorf.fna", "w") as fn, open(p + ".smorf.gff", "w") as fg, \
            open(p + ".smorf_table.tsv", "w") as ft:
        fg.write("##gff-version 3\n")
        ft.write("\t".join(cols) + "\n")
        for x in smorfs:
            h = "%s %s:%d-%d(%s) set=%s start=%s rbs=%s" % (x["id"], x["contig"], x["start"], x["end"], x["strand"],
                                                          x["set"], x["start_status"], x["rbs_status"])
            fa.write(">%s\n%s\n" % (h, x["prot"]))
            fn.write(">%s\n%s\n" % (h, dna_of(seq, x)))
            fg.write("%s\tsmORF_module\tCDS\t%d\t%d\t.\t%s\t0\tID=%s;%s\n"
                     % (x["contig"], x["start"], x["end"], x["strand"], x["id"], gattr(x)))
            row = dict(x, smORF=x["id"], sample=n, length_aa=len(x["prot"]), imp3_locus_tag=x["locus_tag"],
                       mg_reads=("%.2f" % x["mg_reads"]) if "mg_reads" in x else "",
                       mt_reads=("%.2f" % x["mt_reads"]) if "mt_reads" in x else "", support=x.get("support", ""))
            ft.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
    with open(p + ".smorf_stats.tsv", "w") as fh:
        fh.write("\t".join(stats) + "\n" + "\t".join(str(v) for v in stats.values()) + "\n")

    # IMP3 GFF with the smORFs: IMP3 CDS that are smORFs tagged, new smORFs added; in assembly order
    imp3 = {(x["contig"], x["start"], x["end"], x["strand"]): x for x in smorfs if x["set"] == "imp3"}
    new = collections.defaultdict(list)
    for x in smorfs:
        if x["set"] != "imp3":
            new[x["contig"]].append([x["contig"], "smORF_module", "CDS", x["start"], x["end"], ".", x["strand"], "0",
                                     "ID=%s;product=small protein (predicted smORF);%s" % (x["id"], gattr(x))])
    rows = collections.defaultdict(list)
    for r in gff_rows:
        r = list(r)
        x = imp3.get((r[0], r[3], r[4], r[6])) if r[2] == "CDS" else None
        if x is not None:
            r[8] = r[8].rstrip(";") + ";smORF=true;" + gattr(x)
        rows[r[0]].append(r)
    tmp = dst + ".tmp.gz"
    with gzip.open(tmp, "wt") as fh:
        for h in gff_head:
            fh.write(h + "\n")
        for c in list(seq) + [c for c in rows if c not in seq]:
            for r in sorted(rows.get(c, []) + new.get(c, []), key=lambda r: (r[3], r[4])):
                fh.write("\t".join(str(v) for v in r) + "\n")
    os.replace(tmp, dst)
    log("written: %s.smorf.{faa,fna,gff}, smorf_table.tsv, smorf_stats.tsv; %s" % (p, dst))


if __name__ == "__main__":
    main()
