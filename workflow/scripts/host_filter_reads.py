#!/usr/bin/env python3
"""Keep the reads that do not align to a filter reference (PhiX, extra filter libraries), from bwa-mem2 SAM on stdin
(primary records only: samtools view -F 2304; mates adjacent as bwa-mem2 writes them).

A read counts as aligned only when >= --min-aligned-fraction of its length aligns (M/=/X bases), >= --min-aligned-length
bases, and the mismatches (NM) are <= --max-mismatch-fraction of the aligned bases. Paired mode removes a pair when
either mate counts and writes the other pairs to --out1 / --out2 (no singletons are lost); single-end mode writes the
reads that do not count to --out. Reverse-strand records are written back in their sequenced orientation.

bwa-mem2 reports a read as mapped for any local alignment above its minimum score; with short-insert libraries its
pairing and mate rescue then mark many microbial pairs as mapped (Expobiome: 12-31 % of the pairs against T2T-CHM13,
of which < 6 % human). The length-relative rule removes these (0.01-0.1 % left, as Deacon and Hostile).
"""
import argparse, gzip, json, re, sys

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse(line, a):
    f = line.rstrip("\n").split("\t")
    flag = int(f[1])
    seq, qual = f[9], f[10]
    if flag & 16:
        seq, qual = seq.translate(RC)[::-1], qual[::-1]
    hit = False
    if not flag & 4:
        al = sum(int(n) for n, op in CIG.findall(f[5]) if op in "M=X")
        nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 0)
        hit = al >= a.min_aligned_length and al >= a.min_aligned_fraction * len(seq) and nm <= a.max_mismatch_fraction * al
    return f[0], hit, seq, qual, bool(flag & 64), bool(flag & 128)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paired", action="store_true")
    p.add_argument("--out1"); p.add_argument("--out2"); p.add_argument("--out")
    p.add_argument("--stats")
    p.add_argument("--min-aligned-fraction", type=float, default=0.8)
    p.add_argument("--min-aligned-length", type=int, default=30)
    p.add_argument("--max-mismatch-fraction", type=float, default=0.05)
    a = p.parse_args()
    n = {"records": 0, "pairs" if a.paired else "reads": 0, "removed": 0, "rule": {
        "min_aligned_fraction": a.min_aligned_fraction, "min_aligned_length": a.min_aligned_length,
        "max_mismatch_fraction": a.max_mismatch_fraction}}
    if a.paired:
        o1, o2 = gzip.open(a.out1, "wt", compresslevel=1), gzip.open(a.out2, "wt", compresslevel=1)
        buf = None
        for line in sys.stdin:
            if line.startswith("@"):
                continue
            r = parse(line, a); n["records"] += 1
            if buf is None or buf[0] != r[0]:
                if buf is not None:
                    raise SystemExit("unpaired record %s: input must keep mates adjacent" % buf[0])
                buf = r
                continue
            r1, r2 = (buf, r) if buf[4] or r[5] else (r, buf)
            buf = None
            n["pairs"] += 1
            if r1[1] or r2[1]:
                n["removed"] += 1
                continue
            o1.write("@%s/1\n%s\n+\n%s\n" % (r1[0], r1[2], r1[3]))
            o2.write("@%s/2\n%s\n+\n%s\n" % (r2[0], r2[2], r2[3]))
        if buf is not None:
            raise SystemExit("unpaired record %s at the end of the input" % buf[0])
        o1.close(); o2.close()
    else:
        o = gzip.open(a.out, "wt", compresslevel=1)
        for line in sys.stdin:
            if line.startswith("@"):
                continue
            r = parse(line, a); n["records"] += 1; n["reads"] += 1
            if r[1]:
                n["removed"] += 1
                continue
            o.write("@%s\n%s\n+\n%s\n" % (r[0], r[2], r[3]))
        o.close()
    tot = n["pairs"] if a.paired else n["reads"]
    n["removed_fraction"] = n["removed"] / tot if tot else 0.0
    if a.stats:
        json.dump(n, open(a.stats, "w"), indent=1)
    sys.stderr.write("host_filter_reads: %d of %d %s removed (%.3f %%)\n" % (n["removed"], tot, "pairs" if a.paired
                     else "reads", 100 * n["removed_fraction"]))


if __name__ == "__main__":
    main()
