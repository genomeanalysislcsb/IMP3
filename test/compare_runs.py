#!/usr/bin/env python3
"""Compare two IMP3 output folders of the same sample, e.g. a native run and a container run.

  python test/compare_runs.py REFERENCE_DIR NEW_DIR [--top 15] [--out report.tsv]

Compares: status markers, files (missing / new, size changes of large outputs), assembly (contigs, length, N50),
bins per binner (number, total size, size distribution),
Stats/all_stats.tsv, final bins (CheckM2 quality tiers), Bracken profiles per rank (Bray-Curtis, top taxa),
Mantis featureCounts per database (terms, reads, Spearman). Needs Python 3 with pandas.
Small differences are expected (multi-threaded assemblers and binners are not deterministic); look for missing
steps or files and for large changes."""
import argparse
import glob
import gzip
import os
import re
import sys

import pandas as pd

SKIP_DIRS = (".snakemake", "logs", "logs_cluster_jobs", "config", "benchmarks")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def section(t):
    print("\n== %s" % t)


def files(d):
    out = {}
    for root, dirs, fs in os.walk(d):
        dirs[:] = [x for x in dirs if x not in SKIP_DIRS and not x.startswith(".")]
        for f in fs:
            p = os.path.join(root, f)
            try:
                out[os.path.relpath(p, d)] = os.path.getsize(p)
            except OSError:
                pass
    return out


def fasta_stats(p):
    lens, n = [], 0
    with opener(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    lens.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n:
        lens.append(n)
    lens.sort(reverse=True)
    tot, acc, n50 = sum(lens), 0, 0
    for x in lens:
        acc += x
        if acc >= tot / 2:
            n50 = x
            break
    return {"contigs": len(lens), "total_bp": tot, "N50": n50, "max": lens[0] if lens else 0,
            "contigs_>=1kb": sum(1 for x in lens if x >= 1000)}


def first(d, patterns):
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(d, pat)) + glob.glob(os.path.join(d, pat + ".gz")))
        if hits:
            return hits[0]
    return None


def rel(a, b):
    try:
        return "%+.1f%%" % (100.0 * (b - a) / a) if a else ("=" if a == b else "new")
    except TypeError:
        return ""


def bray_curtis(a, b):
    idx = a.index.union(b.index)
    a, b = a.reindex(idx, fill_value=0), b.reindex(idx, fill_value=0)
    s = (a + b).sum()
    return float((a - b).abs().sum() / s) if s else 0.0


def tiers(df):
    comp = next((c for c in df.columns if "ompleteness" in c), None)
    cont = next((c for c in df.columns if "ontamination" in c), None)
    if comp is None or cont is None:
        return None
    c, k = pd.to_numeric(df[comp], errors="coerce"), pd.to_numeric(df[cont], errors="coerce")
    return {"bins": len(df), "HQ (>90/<5)": int(((c > 90) & (k < 5)).sum()),
            "MQ (>=50/<10)": int(((c >= 50) & (k < 10) & ~((c > 90) & (k < 5))).sum()),
            "LQ": int((~((c >= 50) & (k < 10))).sum())}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("--top", type=int, default=15, help="top taxa per profile to show")
    ap.add_argument("--min-size", type=float, default=1.0, help="MB: size changes reported for files >= this")
    ap.add_argument("--size-change", type=float, default=0.2, help="report size changes larger than this fraction")
    a = ap.parse_args(argv)
    R, N = os.path.abspath(a.ref), os.path.abspath(a.new)
    print("reference: %s\nnew:       %s" % (R, N))
    problems = []

    section("status markers")
    sr = {os.path.basename(p) for p in glob.glob(os.path.join(R, "status", "*"))}
    sn = {os.path.basename(p) for p in glob.glob(os.path.join(N, "status", "*"))}
    print("both: %s" % (", ".join(sorted(sr & sn)) or "-"))
    if sr - sn:
        print("ONLY in reference: %s" % ", ".join(sorted(sr - sn)))
        problems.append("status markers missing in new: %s" % ", ".join(sorted(sr - sn)))
    if sn - sr:
        print("only in new: %s" % ", ".join(sorted(sn - sr)))

    section("files (without %s)" % ", ".join(SKIP_DIRS))
    fr, fn = files(R), files(N)

    def is_bin(p):          # bin FASTAs: names and sizes differ between runs (binners are not deterministic)
        return p.startswith("Binning/") and "/bins/" in p
    print("%d files in reference, %d in new (bin FASTAs: per binner, see below)" % (len(fr), len(fn)))
    miss = sorted(p for p in set(fr) - set(fn) if not is_bin(p))
    new = sorted(p for p in set(fn) - set(fr) if not is_bin(p))
    for p in miss[:60]:
        print("  missing  %s  (%.1f MB)" % (p, fr[p] / 1e6))
    if len(miss) > 60:
        print("  ... %d more" % (len(miss) - 60))
    if miss:
        problems.append("%d files missing in new: %s" % (len(miss), ", ".join(miss[:5]) + (" ..." if len(miss) > 5 else "")))
    for p in new[:30]:
        print("  new      %s" % p)
    if len(new) > 30:
        print("  ... %d more new" % (len(new) - 30))
    binners = sorted({p.split("/")[1] for p in list(fr) + list(fn) if is_bin(p)})
    bin_rows = []
    for b in binners:
        row = [b]
        for d, f in ((R, fr), (N, fn)):
            sizes = [fasta_stats(os.path.join(d, p))["total_bp"] for p in sorted(f)
                     if is_bin(p) and p.split("/")[1] == b and re.search(r"\.(fa|fasta|fna)(\.gz)?$", p)
                     and p.count("/") == 3]          # Binning/<binner>/bins/<bin>.fa(.gz), not per-bin subfolders
            row.append(sizes)
        bin_rows.append(row)
        if row[1] and not row[2]:
            problems.append("no %s bins in new" % b)
    big = [(p, fr[p], fn[p]) for p in set(fr) & set(fn) if not is_bin(p) and max(fr[p], fn[p]) >= a.min_size * 1e6
           and abs(fn[p] - fr[p]) > a.size_change * max(fr[p], 1)]
    for p, x, y in sorted(big):
        print("  size     %s  %.1f -> %.1f MB (%s)" % (p, x / 1e6, y / 1e6, rel(x, y)))

    section("assembly")
    pr = first(R, ["Assembly/*.assembly.merged.fa"])
    pn = first(N, ["Assembly/*.assembly.merged.fa"])
    if pr and pn:
        x, y = fasta_stats(pr), fasta_stats(pn)
        for k in x:
            print("  %-14s %12s %12s  %s" % (k, x[k], y[k], rel(x[k], y[k])))
    else:
        print("  assembly missing: %s" % ("reference" if not pr else "new"))
        if not pn:
            problems.append("no assembly in new")

    section("Stats/all_stats.tsv")
    ar, an = first(R, ["Stats/all_stats.tsv"]) or "-", first(N, ["Stats/all_stats.tsv"]) or "-"
    if os.path.exists(ar) and os.path.exists(an):
        x, y = pd.read_csv(ar, sep="\t"), pd.read_csv(an, sep="\t")
        keys = [c for c in ("type", "step") if c in x.columns and c in y.columns]
        val = "number" if "number" in x.columns else x.columns[-1]
        m = x.merge(y, on=keys, how="outer", suffixes=("_ref", "_new"))
        m["change"] = [rel(p, q) if pd.notna(p) and pd.notna(q) else "missing" for p, q in
                       zip(m[val + "_ref"], m[val + "_new"])]
        with pd.option_context("display.max_rows", 500, "display.width", 200):
            print(m[keys + [val + "_ref", val + "_new", "change"]].to_string(index=False))
    else:
        print("  not in %s" % ", ".join(n for n, p in (("reference", ar), ("new", an)) if not os.path.exists(p)))

    section("bins per binner (bin names and contents differ between runs: number and size distribution)")
    if bin_rows:
        print("  %-12s %17s %17s   %-41s %17s" % ("", "bins", "total Mb", "median Mb (IQR)", "bins >= 1 Mb"))
        print("  %-12s %8s %8s %8s %8s   %-20s %-20s %8s %8s" % ("binner", "ref", "new", "ref", "new", "ref", "new",
                                                              "ref", "new"))
        for b, x, y in bin_rows:
            def q(v, f):
                return pd.Series(v).quantile(f) / 1e6 if v else 0
            print("  %-12s %8d %8d %8.1f %8.1f   %-20s %-20s %8d %8d" % (
                b, len(x), len(y), sum(x) / 1e6, sum(y) / 1e6,
                "%.2f (%.2f-%.2f)" % (q(x, .5), q(x, .25), q(x, .75)), "%.2f (%.2f-%.2f)" % (q(y, .5), q(y, .25), q(y, .75)),
                sum(1 for v in x if v >= 1e6), sum(1 for v in y if v >= 1e6)))
    else:
        print("  no Binning/<binner>/bins/ folders")

    section("final bins (CheckM2 tiers)")
    pats = ["Stats/all_bin_stats.tsv", "Stats/*/*.bins.tsv", "Binning/MetaWrap/*checkm2*/quality_report.tsv",
            "Binning/*/*checkm2*/quality_report.tsv", "Binning/*/quality_report.tsv"]
    br, bn = first(R, pats), first(N, pats)
    t = {}
    for name, p in (("reference", br), ("new", bn)):
        if p:
            try:
                t[name] = tiers(pd.read_csv(p, sep="\t"))
            except Exception as e:
                t[name] = None
                print("  %s: cannot read %s (%s)" % (name, p, e))
            print("  %-9s %s" % (name, os.path.relpath(p, R if name == "reference" else N)))
        else:
            print("  %-9s no bin table (%s)" % (name, ", ".join(pats)))
    if t.get("reference") or t.get("new"):
        keys = list((t.get("reference") or t.get("new")).keys())
        for k in keys:
            print("  %-14s %6s %6s" % (k, (t.get("reference") or {}).get(k, "-"), (t.get("new") or {}).get(k, "-")))

    section("Bracken profiles (fraction of reads)")
    if not glob.glob(os.path.join(R, "Analysis", "taxonomy", "bracken", "*", "*.bracken*")):
        print("  no Analysis/taxonomy/bracken/<omic>/<rank>.bracken in the reference; found: %s" % (", ".join(
            os.path.relpath(p, R) for p in sorted(glob.glob(os.path.join(R, "Analysis", "taxonomy", "*", "*bracken*")))[:8])
            or "-"))
    for pr in sorted(glob.glob(os.path.join(R, "Analysis", "taxonomy", "bracken", "*", "*.bracken*"))):
        r = os.path.relpath(pr, R)
        pn = first(N, [r[:-3] if r.endswith(".gz") else r]) or os.path.join(N, r)
        if not os.path.exists(pn):
            print("  %s: missing in new" % r)
            continue
        try:
            x = pd.read_csv(pr, sep="\t").set_index("name")["fraction_total_reads"]
            y = pd.read_csv(pn, sep="\t").set_index("name")["fraction_total_reads"]
        except Exception as e:
            print("  %s: cannot read (%s)" % (r, e))
            continue
        print("  %-40s taxa %5d %5d  Bray-Curtis %.3f" % (r, len(x), len(y), bray_curtis(x, y)))
        if r.endswith(("/G.bracken", "/S.bracken", "/G.bracken.gz", "/S.bracken.gz")):
            t = pd.DataFrame({"ref": x, "new": y}).fillna(0).sort_values("ref", ascending=False).head(a.top)
            for name, row in t.iterrows():
                print("      %-40s %7.4f %7.4f" % (str(name)[:40], row["ref"], row["new"]))

    section("Mantis featureCounts per database")
    if not glob.glob(os.path.join(R, "Analysis", "annotation", "*.mantis.*.counts.tsv*")):
        print("  no Analysis/annotation/<omic>.mantis.<db>.counts.tsv in the reference; found: %s" % (", ".join(
            os.path.relpath(p, R) for p in sorted(glob.glob(os.path.join(R, "Analysis", "annotation", "*counts*")))[:8])
            or "-"))
    for pr in sorted(glob.glob(os.path.join(R, "Analysis", "annotation", "*.mantis.*.counts.tsv*"))):
        r = os.path.relpath(pr, R)
        pn = first(N, [r[:-3] if r.endswith(".gz") else r]) or os.path.join(N, r)
        if not os.path.exists(pn):
            print("  %s: missing in new" % r)
            continue
        try:
            x = pd.to_numeric(pd.read_csv(pr, sep="\t", comment="#").set_index("Geneid").iloc[:, -1], errors="coerce")
            y = pd.to_numeric(pd.read_csv(pn, sep="\t", comment="#").set_index("Geneid").iloc[:, -1], errors="coerce")
        except Exception as e:
            print("  %s: cannot read (%s)" % (r, e))
            continue
        both = pd.DataFrame({"ref": x, "new": y}).fillna(0)
        rho = both["ref"].corr(both["new"], method="spearman") if len(both) > 2 else float("nan")
        print("  %-40s terms %6d %6d  reads %10d %10d (%s)  Spearman %.3f" % (
            os.path.basename(r), int((x > 0).sum()), int((y > 0).sum()), int(x.sum()), int(y.sum()),
            rel(x.sum(), y.sum()), rho))

    section("summary")
    print("\n".join("  - " + p for p in problems) if problems else "  no missing steps or files")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
