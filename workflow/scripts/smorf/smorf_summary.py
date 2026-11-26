#!/usr/bin/env python3
"""Numbers of the smORF module of IMP3.

  sample  one sample: intergenic regions and IMP3 CDS <= 100 aa (regions: intergenic), genes of the two pyrodigal-gv
          runs, smORF candidates, after deduplication, hits of each
          false-positive filter (Rfam, AntiFam, Aragorn), final smORFs and their lengths
  cohort  all samples: the per-sample table, and per clustering level the clusters with their samples
          (<prefix>_rep_samples.tsv) and the numbers of the paper (clusters, singletons, unique to one sample,
          in >= 2 samples, core = in every sample)
"""
import argparse
import statistics
from collections import defaultdict


def fasta_lengths(path):
    out, name, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip().replace("*", ""))
    if name is not None:
        out[name] = n
    return out


def count_headers(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def rfam_ids(path):
    s = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                p = line.split()
                if len(p) > 3:
                    s.add(p[3])
    return s


def antifam_ids(path):
    s = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                p = line.split()
                if len(p) > 3:
                    s.add(p[0])
    return s


def aragorn_ids(path):
    s, lines = set(), open(path).read().splitlines()
    for i, line in enumerate(lines):
        if line.startswith(">") and "::" in line and i + 1 < len(lines) and "0 genes found" not in lines[i + 1]:
            s.add(line[1:].split(" ")[0])
    return s


def cmd_sample(a):
    final = fasta_lengths(a.final)
    rf, af = rfam_ids(a.rfam), antifam_ids(a.antifam)
    ar = aragorn_ids(a.aragorn_std) | aragorn_ids(a.aragorn_bact)
    lens = list(final.values())
    regions = {}
    if a.regions:
        h, v = open(a.regions).read().splitlines()[:2]
        regions = dict(zip(h.split("\t"), v.split("\t")))
    imp3 = count_headers(a.cand_imp3) if a.cand_imp3 else 0
    n_imp3 = sum(1 for k in final if "::imp3_" in k)
    row = [("sample", a.sample),
           ("intergenic_regions", regions.get("intergenic_regions", "NA")),
           ("intergenic_bp", regions.get("intergenic_bp", "NA")),
           ("assembly_bp", regions.get("assembly_bp", "NA")),
           ("smorf_candidates_imp3", imp3 if a.cand_imp3 else "NA"),
           ("genes_default", count_headers(a.genes_default)), ("genes_small", count_headers(a.genes_small)),
           ("smorf_candidates_default", count_headers(a.cand_default)),
           ("smorf_candidates_small", count_headers(a.cand_small)),
           ("smorfs_deduplicated", count_headers(a.potential)),
           ("removed_rfam", len(rf)), ("removed_antifam", len(af)), ("removed_aragorn", len(ar)),
           ("removed_total", len(rf | af | ar)),
           ("smorfs", len(final)),
           ("smorfs_imp3_cds", n_imp3), ("smorfs_new", len(final) - n_imp3),
           ("smorfs_le50aa", sum(1 for x in lens if x <= 50)),
           ("smorfs_lt20aa", sum(1 for x in lens if x < 20)),
           ("length_min", min(lens) if lens else 0), ("length_median", statistics.median(lens) if lens else 0),
           ("length_max", max(lens) if lens else 0)]
    with open(a.out, "w") as fo:
        fo.write("\t".join(k for k, _ in row) + "\n")
        fo.write("\t".join(str(v) for _, v in row) + "\n")


def cmd_cohort(a):
    # per-sample rows
    rows = []
    for f in a.samples:
        with open(f) as fh:
            lines = fh.read().splitlines()
        hdr, rows_ = lines[0], lines[1:]
        rows += rows_
    with open(a.out, "w") as fo:
        fo.write(hdr + "\n" + "\n".join(rows) + "\n")
    n_samples = len(rows)
    # clustering levels
    out = []
    for lvl in a.clusters:
        label, tsv, rep = lvl.split(",")
        lengths = fasta_lengths(rep)
        members = defaultdict(list)
        with open(tsv) as fh:
            for line in fh:
                r, m = line.rstrip("\n").split("\t")
                members[r].append(m)
        n_single = n_unique = n_multi_sample = n_core = 0
        with open(rep.replace("_rep_seq.fasta", "_rep_samples.tsv"), "w") as fo:
            fo.write("representative\tlength_aa\tmembers\tsamples\tsample_names\n")
            for r, ms in members.items():
                smp = sorted({m.split("::", 1)[0] for m in ms})
                n_single += len(ms) == 1
                n_unique += len(smp) == 1
                n_multi_sample += len(smp) >= 2
                n_core += len(smp) == n_samples
                fo.write("%s\t%d\t%d\t%d\t%s\n" % (r, lengths.get(r, 0), len(ms), len(smp), ",".join(smp)))
        n_all = sum(len(v) for v in members.values())
        out.append((label, n_all, len(members), n_single, n_unique, n_multi_sample, n_core))
    with open(a.cluster_summary, "w") as fo:
        fo.write("level\tsmORFs\tclusters\tsingletons\tin_one_sample\tin_2plus_samples\tcore_all_%d_samples\n"
                 % n_samples)
        for r in out:
            fo.write("\t".join(str(x) for x in r) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("sample")
    for k in ("sample", "genes_default", "genes_small", "cand_default", "cand_small", "potential", "rfam", "antifam",
              "aragorn_std", "aragorn_bact", "final", "out"):
        p.add_argument("--" + k.replace("_", "-"), dest=k, required=True)
    p.add_argument("--cand-imp3", default="", help="IMP3 CDS <= max_aa (regions: intergenic)")
    p.add_argument("--regions", default="", help="Stats/<sample>.intergenic.tsv (regions: intergenic)")
    p = sub.add_parser("cohort")
    p.add_argument("--samples", nargs="+", required=True, help="per-sample tables")
    p.add_argument("--clusters", nargs="*", default=[], help="label,cluster.tsv,rep_seq.fasta")
    p.add_argument("--out", required=True)
    p.add_argument("--cluster-summary", required=True)
    a = ap.parse_args()
    cmd_sample(a) if a.cmd == "sample" else cmd_cohort(a)


if __name__ == "__main__":
    main()
