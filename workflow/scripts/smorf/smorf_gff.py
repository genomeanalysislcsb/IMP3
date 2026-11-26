#!/usr/bin/env python3
"""IMP3 annotation with the smORFs, and one table of all smORFs of a sample (smORF module of IMP3).

GFF (--out; a new file, the IMP3 output stays unchanged):
  - an IMP3 CDS that is a final smORF (same contig, strand and stop codon) gets smORF=true; smorf_id=...
  - every other smORF is added as a CDS line, source smORF_module, ID / locus_tag <sample>_smORF_<n>,
    product=small protein (predicted smORF)
  - smORF attributes: smorf_status=predicted|validated (metaproteomics), smorf_length_aa, smorf_set
    (imp3 / default / small), smorf_overlaps (IMP3 features other than the CDS itself that overlap the smORF, e.g.
    regulatory_region:<product> for leader peptides in riboswitches / attenuators, ncRNA:<product>)
Table (--table): smORF, sample, set, contig, start, end, strand, length_aa, imp3_locus_tag, gff_id, status,
unique_peptides, overlaps, and the read columns of the read support (mg_reads, mt_reads, support) when given.
"""
import argparse
import csv
import gzip
from collections import defaultdict


def zopen(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def stop_key(contig, start, end, strand):
    return (contig, strand, int(end) if strand == "+" else int(start))


def attr(a, key):
    return next((x.split("=", 1)[1] for x in a.split(";") if x.startswith(key + "=")), "")


def read_table(path, key):
    if not path:
        return {}
    with open(path) as fh:
        return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imp3-gff", required=True)
    ap.add_argument("--smorf-gff", required=True, help="smORF/<sample>.smorf.gff of the module")
    ap.add_argument("--smorf-faa", required=True, help="smORF/<sample>.smorf.faa (lengths)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--status", default="", help="metaproteomics validation table (smorf_metap.py)")
    ap.add_argument("--reads", default="", help="read support table (smorf_read_support.py summary)")
    ap.add_argument("--overlap-types", default="regulatory_region ncRNA tRNA tmRNA rRNA CRISPR CDS")
    ap.add_argument("--out", required=True)
    ap.add_argument("--table", required=True)
    a = ap.parse_args()

    length, name = {}, None
    for line in open(a.smorf_faa):
        if line.startswith(">"):
            name = line[1:].split()[0]
            length[name] = 0
        else:
            length[name] += len(line.strip())
    status = read_table(a.status, "smORF")
    reads = read_table(a.reads, "smORF")

    smorfs = []                                  # (contig, start, end, strand, id, set)
    for line in open(a.smorf_gff):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) != 9:
            continue
        sset, contig = p[0].split("_", 1)
        smorfs.append((contig, int(p[3]), int(p[4]), p[6], attr(p[8], "ID"), sset))

    header, lines, contig_order = [], [], {}
    feats = defaultdict(list)                    # contig -> [(start, end, type, label, stop_key)]
    otypes = set(a.overlap_types.split())
    with zopen(a.imp3_gff) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                if not lines:
                    header.append(line)
                continue
            p = line.split("\t")
            if len(p) != 9:
                continue
            contig_order.setdefault(p[0], len(contig_order))
            lines.append(p)
            if p[2] in otypes:
                label = attr(p[8], "product") or attr(p[8], "Name") or attr(p[8], "ID")
                feats[p[0]].append((int(p[3]), int(p[4]), p[2], label.replace(";", ","),
                                    stop_key(p[0], p[3], p[4], p[6]) if p[2] == "CDS" else None,
                                    attr(p[8], "locus_tag") or attr(p[8], "ID")))

    def overlaps(c, s, e, own):
        out = []
        for fs, fe, t, label, k, _ in feats.get(c, ()):
            if fe < s or fs > e or (k is not None and k == own):
                continue
            out.append("%s:%s" % (t, label.replace(" ", "_").replace(",", "_")))
        return ",".join(sorted(set(out)))

    # IMP3 CDS that are smORFs: same stop codon
    cds_by_stop = {}
    for c, lst in feats.items():
        for fs, fe, t, label, k, tag in lst:
            if t == "CDS":
                cds_by_stop[k] = tag
    rows, imp3_tag = {}, {}
    for c, s, e, st, sid, sset in smorfs:
        k = stop_key(c, s, e, st)
        rows[sid] = dict(contig=c, start=s, end=e, strand=st, set=sset, key=k, overlaps=overlaps(c, s, e, k),
                         imp3=cds_by_stop.get(k, ""))
    first_on_cds = {}
    for sid, r in rows.items():
        if r["imp3"]:
            first_on_cds.setdefault(r["key"], sid)

    def smorf_attrs(sid):
        r, stt = rows[sid], status.get(sid, {})
        x = "smORF=true;smorf_id=%s;smorf_set=%s;smorf_length_aa=%d;smorf_status=%s" % (
            sid, r["set"], length.get(sid, 0), stt.get("status", "predicted"))
        if stt.get("unique_peptides", "0") not in ("", "0"):
            x += ";smorf_unique_peptides=%s" % stt["unique_peptides"]
        if r["overlaps"]:
            x += ";smorf_overlaps=%s" % r["overlaps"]
        return x

    gff_id = {}
    for p in lines:
        if p[2] == "CDS":
            k = stop_key(p[0], p[3], p[4], p[6])
            if k in first_on_cds:
                sid = first_on_cds[k]
                p[8] = p[8].rstrip(";") + ";" + smorf_attrs(sid)
                gff_id[sid] = attr(p[8], "ID")
    n_new = 0
    for c, s, e, st, sid, sset in sorted(smorfs, key=lambda r: (r[0], r[1])):
        if sid in gff_id:
            continue
        n_new += 1
        contig_order.setdefault(c, len(contig_order))
        tag = "%s_smORF_%d" % (a.sample, n_new)
        gff_id[sid] = tag
        lines.append([c, "smORF_module", "CDS", str(s), str(e), ".", st, "0",
                      "ID=%s;locus_tag=%s;product=small protein (predicted smORF);%s" % (tag, tag, smorf_attrs(sid))])
    lines.sort(key=lambda p: (contig_order[p[0]], int(p[3]), int(p[4])))
    with open(a.out, "w") as fo:
        fo.write("\n".join(header or ["##gff-version 3"]) + "\n")
        for p in lines:
            fo.write("\t".join(p) + "\n")

    rcols = [c for c in ("mg_reads", "mt_reads", "support") if reads and c in next(iter(reads.values()))]
    with open(a.table, "w") as fo:
        fo.write("\t".join(["smORF", "sample", "set", "contig", "start", "end", "strand", "length_aa",
                            "imp3_locus_tag", "gff_id", "status", "unique_peptides", "overlaps"] + rcols) + "\n")
        for c, s, e, st, sid, sset in sorted(smorfs, key=lambda r: (contig_order.get(r[0], 0), r[1])):
            r, stt, rd = rows[sid], status.get(sid, {}), reads.get(sid, {})
            fo.write("\t".join([sid, a.sample, sset, c, str(s), str(e), st, str(length.get(sid, 0)), r["imp3"],
                                gff_id[sid], stt.get("status", "predicted"), stt.get("unique_peptides", "0"),
                                r["overlaps"]] + [rd.get(x, "") for x in rcols]) + "\n")
    n_val = sum(1 for sid in rows if status.get(sid, {}).get("status") == "validated")
    print("%s\t%d smORFs\t%d IMP3 CDS tagged\t%d added\t%d validated" % (
        a.sample, len(smorfs), len(smorfs) - n_new, n_new, n_val))


if __name__ == "__main__":
    main()
