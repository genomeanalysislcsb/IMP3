#!/usr/bin/env python3
"""Remapped reads per MAG (bins) of another IMP3 run.

  remap_bins.py abundance --coverage <samtools coverage> --bins <contig2bin> --out <tsv>
      per bin: contigs, length, reads, share of the mapped reads, length-weighted mean depth (unbinned contigs: bin
      'unbinned')
  remap_bins.py functions --counts <featureCounts -f per-contig table> --bins <contig2bin> --out <tsv>
      per bin and term (featureCounts Geneid, e.g. a KO): counts summed over the bin's contigs
"""
import argparse, collections, gzip, sys


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def contig2bin(path):
    c2b = {}
    for line in opener(path):
        x = line.split()
        if len(x) >= 2 and not line.startswith("#"):
            c2b[x[0]] = x[1]
    return c2b


def abundance(a):
    c2b = contig2bin(a.bins)
    agg = collections.defaultdict(lambda: [0, 0, 0, 0.0])          # contigs, length, reads, depth x length
    total = 0
    for line in opener(a.coverage):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")              # rname startpos endpos numreads covbases coverage meandepth ...
        ln, reads, depth = int(f[2]) - int(f[1]) + 1, int(f[3]), float(f[6])
        b = agg[c2b.get(f[0], "unbinned")]
        b[0] += 1; b[1] += ln; b[2] += reads; b[3] += depth * ln
        total += reads
    with open(a.out, "w") as o:
        o.write("bin\tcontigs\tlength\treads\tshare_mapped\tmean_depth\n")
        for b, (n, ln, r, dl) in sorted(agg.items(), key=lambda t: -t[1][2]):
            o.write("%s\t%d\t%d\t%d\t%.6f\t%.4f\n" % (b, n, ln, r, r / total if total else 0, dl / ln if ln else 0))


def functions(a):
    c2b = contig2bin(a.bins)
    agg = collections.defaultdict(float)
    head = None
    for line in opener(a.counts):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if head is None:
            head = f                                    # Geneid Chr Start End Strand Length <bam>
            continue
        agg[(c2b.get(f[1], "unbinned"), f[0])] += float(f[6])
    with open(a.out, "w") as o:
        o.write("bin\tterm\tcount\n")
        for (b, t), c in sorted(agg.items()):
            if c > 0:
                o.write("%s\t%s\t%.4f\n" % (b, t, c))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest="cmd", required=True)
    x = s.add_parser("abundance"); x.add_argument("--coverage", required=True)
    y = s.add_parser("functions"); y.add_argument("--counts", required=True)
    for q in (x, y):
        q.add_argument("--bins", required=True); q.add_argument("--out", required=True)
    a = p.parse_args()
    (abundance if a.cmd == "abundance" else functions)(a)


if __name__ == "__main__":
    main()
