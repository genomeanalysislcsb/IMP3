#!/usr/bin/env python3
"""smORF candidates from the IMP3 annotation (smORF module of IMP3, prediction: regions: intergenic).

Instead of calling genes again on the whole assembly (the paper's pipeline), the module uses the gene calls of the
IMP3 run and looks for new smORFs only where IMP3 has no feature:

  regions  intergenic regions of the IMP3 annotation (GFF; feature types in --occupied, default CDS, tRNA, tmRNA,
           rRNA, ncRNA, CRISPR) of at least --min-length nt, each extended by --flank nt into its neighbours; written
           as FASTA fragments for pyrodigal-gv plus a table fragment -> contig, offset
  imp3     the CDS of IMP3 with proteins of <= --max-aa aa (complete genes; not truncated, not at a contig end) as
           smORF candidates in the format of the prediction (pyrodigal_gv/imp3/<sample>.{smorf.faa,fna,gff})
  lift     pyrodigal-gv output on the fragments -> contig coordinates, genes named <contig>_<k> (k by position), so
           that the rest of the pipeline works as on whole contigs; genes extending beyond the fragment core into a
           neighbouring IMP3 feature by more than the flank are impossible by construction
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def zopen(path, mode="rt"):
    if "w" in mode:
        return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, mode) if gz else open(path, mode)


def read_fasta(path):
    name, seq = None, []
    with zopen(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def attrs(s):
    out = {}
    for x in s.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def read_gff(path):
    """[(contig, type, start, end, strand, attributes dict)] (1-based, inclusive)"""
    out = []
    with zopen(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                continue
            out.append((p[0], p[2], int(p[3]), int(p[4]), p[6], attrs(p[8])))
    return out


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def cmd_regions(a):
    contigs = {n.split()[0]: s for n, s in read_fasta(a.assembly)}
    occupied = set(a.occupied.split())
    occ = defaultdict(list)
    for c, t, s, e, _, _ in read_gff(a.gff):
        if t in occupied and c in contigs:
            occ[c].append((s, e))
    n = n_bp = 0
    with zopen(a.out_fasta, "wt") as fo, open(a.out_map, "w") as fm:
        fm.write("fragment\tcontig\toffset\tcore_start\tcore_end\tlength\n")
        for c, seq in contigs.items():
            L, pos = len(seq), 1
            gaps = []
            for s, e in merge(occ.get(c, [])) + [[L + 1, L + 1]]:
                if s - pos >= a.min_length:
                    gaps.append((pos, s - 1))
                pos = max(pos, e + 1)
            for s, e in gaps:
                fs, fe = max(1, s - a.flank), min(L, e + a.flank)
                n += 1
                n_bp += e - s + 1
                frag = "ig%d" % n
                fo.write(">%s\n%s\n" % (frag, seq[fs - 1:fe]))
                fm.write("%s\t%s\t%d\t%d\t%d\t%d\n" % (frag, c, fs - 1, s, e, fe - fs + 1))
    with open(a.stats, "w") as fh:
        fh.write("contigs\tassembly_bp\tintergenic_regions\tintergenic_bp\n%d\t%d\t%d\t%d\n"
                 % (len(contigs), sum(len(s) for s in contigs.values()), n, n_bp))
    sys.stderr.write("%d intergenic regions, %d bp\n" % (n, n_bp))


def cmd_imp3(a):
    contigs = {n.split()[0]: s for n, s in read_fasta(a.assembly)}
    prot = {n.split()[0]: s.replace("*", "") for n, s in read_fasta(a.faa)}
    rows = []
    for c, t, s, e, st, at in read_gff(a.gff):
        if t != "CDS" or c not in contigs:
            continue
        tag = at.get("locus_tag") or at.get("ID")
        p = prot.get(tag) or prot.get(at.get("ID", ""))
        if not p or len(p) > a.max_aa:
            continue
        txt = ";".join("%s=%s" % kv for kv in at.items()).lower()
        if "truncated" in txt or "partial=1" in txt or "partial=01" in txt or "partial=10" in txt:
            continue
        if s <= 3 or e >= len(contigs[c]) - 2:
            continue
        rows.append((c, s, e, st, tag, p))
    rows.sort(key=lambda r: (r[0], r[1]))
    k = defaultdict(int)
    with open(a.out_faa, "w") as fa, open(a.out_fna, "w") as fn, open(a.out_gff, "w") as fg:
        fg.write("##gff-version  3\n")
        for c, s, e, st, tag, p in rows:
            k[c] += 1
            gid = "%s_%d" % (c, k[c])
            dna = contigs[c][s - 1:e]
            if st == "-":
                dna = dna.translate(COMP)[::-1]
            info = "ID=%s;partial=00;locus_tag=%s;source=IMP3" % (gid, tag)
            head = "imp3_%s # %d # %d # %d # %s" % (gid, s, e, 1 if st == "+" else -1, info)
            fa.write(">%s::%s;smORF [len %d]\n%s\n" % (a.sample, head, len(p), p))
            fn.write(">%s\n%s\n" % (head, dna))
            fg.write("imp3_%s\tIMP3\tCDS\t%d\t%d\t.\t%s\t0\t%s;\n" % (c, s, e, st, info))
    sys.stderr.write("%d IMP3 CDS <= %d aa\n" % (len(rows), a.max_aa))


def cmd_lift(a):
    fmap = {}
    with open(a.map) as fh:
        next(fh)
        for line in fh:
            f, c, off, cs, ce, _ = line.rstrip("\n").split("\t")
            fmap[f] = (c, int(off), int(cs), int(ce))
    # genes of the fragments, by position on the contig
    genes = {}
    for c, t, s, e, st, at in read_gff(a.gff_in):
        if c not in fmap:
            continue
        contig, off, cs, ce = fmap[c]
        gs, ge = s + off, e + off
        if ge < cs or gs > ce:          # entirely in the flank: lies in an IMP3 feature
            continue
        genes[at["ID"]] = (contig, gs, ge, st)
    order = sorted(genes, key=lambda g: (genes[g][0], genes[g][1], genes[g][2]))
    k, new = defaultdict(int), {}
    for g in order:
        c = genes[g][0]
        k[c] += 1
        new[g] = "%s_%d" % (c, k[c])
    # GFF
    with open(a.gff_in) as fh, open(a.gff_out, "w") as fo:
        fo.write("##gff-version  3\n")
        lines = {}
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                continue
            g = attrs(p[8]).get("ID")
            if g not in new:
                continue
            c, gs, ge, st = genes[g]
            p[0], p[3], p[4] = c, str(gs), str(ge)
            p[8] = re.sub(r"(^|;)ID=[^;]*", r"\1ID=%s" % new[g], p[8]) + ";fragment=%s" % g.rsplit("_", 1)[0]
            lines[new[g]] = "\t".join(p)
        for g in order:
            fo.write(lines[new[g]] + "\n")
    # FASTA: >ig12_3 # 10 # 120 # 1 # ID=12_3;... -> >contig_k # start # end # strand # ID=contig_k;...
    for fin, fout in ((a.faa_in, a.faa_out), (a.fna_in, a.fna_out)):
        recs = {}
        for h, s in read_fasta(fin):
            g = h.split()[0]
            if g not in new:
                continue
            c, gs, ge, st = genes[g]
            parts = h.split(" # ")
            parts[0] = new[g]
            if len(parts) >= 5:
                parts[1], parts[2] = str(gs), str(ge)
                parts[4] = re.sub(r"^ID=[^;]*", "ID=%s" % new[g], parts[4])
            recs[new[g]] = (" # ".join(parts), s)
        with open(fout, "w") as fo:
            for g in order:
                h, s = recs[new[g]]
                fo.write(">%s\n%s\n" % (h, s))
    sys.stderr.write("%d genes on the intergenic regions\n" % len(order))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("regions")
    p.add_argument("--assembly", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--occupied", default="CDS tRNA tmRNA rRNA ncRNA CRISPR")
    p.add_argument("--min-length", type=int, default=30)
    p.add_argument("--flank", type=int, default=0)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-map", required=True)
    p.add_argument("--stats", required=True)
    p = sub.add_parser("imp3")
    p.add_argument("--sample", required=True)
    p.add_argument("--assembly", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--faa", required=True)
    p.add_argument("--max-aa", type=int, default=100)
    p.add_argument("--out-faa", required=True)
    p.add_argument("--out-fna", required=True)
    p.add_argument("--out-gff", required=True)
    p = sub.add_parser("lift")
    p.add_argument("--map", required=True)
    for k in ("faa", "fna", "gff"):
        p.add_argument("--%s-in" % k, required=True)
        p.add_argument("--%s-out" % k, required=True)
    a = ap.parse_args()
    {"regions": cmd_regions, "imp3": cmd_imp3, "lift": cmd_lift}[a.cmd](a)


if __name__ == "__main__":
    main()
