#!/usr/bin/env python3
"""Metagenomic and metatranscriptomic reads on the predicted smORFs (smORF module of IMP3).

  saf      smORF GFF of the prediction -> SAF for featureCounts on the IMP3 assembly (contig names without the
           default_ / small_ prefix of the two gene-calling runs)
  summary  featureCounts tables (mg, mt) -> one row per smORF: reads per omic and support class
           (MG only / MT only / MG+MT / none), as in the paper (64 % MG only, 11 % MT only, 25 % both)
"""
import argparse
import sys


def saf(gff, out):
    n = 0
    with open(gff) as fh, open(out, "w") as fo:
        fo.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                continue
            chrom = p[0]
            for pre in ("default_", "small_"):
                if chrom.startswith(pre):
                    chrom = chrom[len(pre):]
                    break
            gid = next((x[3:] for x in p[8].split(";") if x.startswith("ID=")), None)
            if gid is None:
                continue
            fo.write("%s\t%s\t%s\t%s\t%s\n" % (gid, chrom, p[3], p[4], p[6]))
            n += 1
    sys.stderr.write("%d smORFs in %s\n" % (n, out))


def read_fc(path):
    """featureCounts table -> {GeneID: (length, reads)}"""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            p = line.rstrip("\n").split("\t")
            out[p[0]] = (int(p[5]), float(p[6]))
    return out


def summary(counts, min_reads, sample, out, stats):
    tabs = {}
    for c in counts:
        omic, path = c.split("=", 1)
        tabs[omic] = read_fc(path)
    omics = [o for o in ("mg", "mt") if o in tabs]
    ids = sorted(set().union(*[set(t) for t in tabs.values()])) if tabs else []
    classes = {"MG only": 0, "MT only": 0, "MG+MT": 0, "none": 0}
    with open(out, "w") as fo:
        fo.write("smORF\tlength_nt\t%s\tsupport\n" % "\t".join("%s_reads" % o for o in omics))
        for i in ids:
            length = next(tabs[o][i][0] for o in omics if i in tabs[o])
            reads = [tabs[o].get(i, (0, 0.0))[1] for o in omics]
            mg = "mg" in omics and tabs["mg"].get(i, (0, 0.0))[1] >= min_reads
            mt = "mt" in omics and tabs["mt"].get(i, (0, 0.0))[1] >= min_reads
            cls = "MG+MT" if mg and mt else "MG only" if mg else "MT only" if mt else "none"
            classes[cls] += 1
            fo.write("%s\t%d\t%s\t%s\n" % (i, length, "\t".join("%.2f" % r for r in reads), cls))
    with open(stats, "w") as fo:
        fo.write("sample\tomics\tsmORFs\t%s\n" % "\t".join(classes))
        fo.write("%s\t%s\t%d\t%s\n" % (sample, "+".join(omics), len(ids), "\t".join(str(v) for v in classes.values())))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("saf")
    p.add_argument("--gff", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("summary")
    p.add_argument("--counts", nargs="+", required=True, help="omic=featureCounts table, e.g. mg=mg.counts.tsv")
    p.add_argument("--min-reads", type=float, default=1)
    p.add_argument("--sample", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--stats", required=True)
    a = ap.parse_args()
    if a.cmd == "saf":
        saf(a.gff, a.out)
    else:
        summary(a.counts, a.min_reads, a.sample, a.out, a.stats)


if __name__ == "__main__":
    main()
