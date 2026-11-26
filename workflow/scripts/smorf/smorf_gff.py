#!/usr/bin/env python3
"""IMP3 annotation with the smORFs (smORF module of IMP3).

Writes the IMP3 GFF of a sample with its final smORFs added, as a new file (the IMP3 output stays unchanged):
  - an IMP3 CDS that is a final smORF (same contig, strand and stop codon) gets smORF=true; smorf_id=<smORF name>
  - every other smORF (new, from the intergenic regions or the gene calling of the module) is added as a CDS line,
    source smORF_module, ID / locus_tag <sample>_smORF_<n>, product=small protein (predicted smORF), with
    smORF=true, smorf_id=<smORF name>, smorf_run=<default|small> and the length in aa
Lines are sorted by contig (order of the IMP3 GFF) and start.
"""
import argparse
import gzip


def zopen(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def stop_key(contig, start, end, strand):
    return (contig, strand, int(end) if strand == "+" else int(start))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imp3-gff", required=True)
    ap.add_argument("--smorf-gff", required=True, help="smORF/<sample>.smorf.gff of the module")
    ap.add_argument("--smorf-faa", required=True, help="smORF/<sample>.smorf.faa (lengths)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    length = {}
    name = None
    for line in open(a.smorf_faa):
        if line.startswith(">"):
            name = line[1:].split()[0]
            length[name] = 0
        else:
            length[name] += len(line.strip())

    smorfs = []                                  # (contig, start, end, strand, id, run)
    for line in open(a.smorf_gff):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) != 9:
            continue
        run, contig = p[0].split("_", 1)
        sid = next((x[3:] for x in p[8].split(";") if x.startswith("ID=")), "")
        smorfs.append((contig, int(p[3]), int(p[4]), p[6], sid, run))

    by_stop = {}                                 # stop codon -> first smORF ending there
    for c, s, e, st, sid, run in smorfs:
        by_stop.setdefault(stop_key(c, s, e, st), (sid, run))
    header, lines, contig_order, matched = [], [], {}, set()
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
            if p[2] == "CDS":
                k = stop_key(p[0], p[3], p[4], p[6])
                if k in by_stop:
                    sid, run = by_stop[k]
                    matched.add(sid)
                    p[8] = p[8].rstrip(";") + ";smORF=true;smorf_id=%s;smorf_length_aa=%d" % (sid, length.get(sid, 0))
            lines.append(p)
    n_new = 0
    for c, s, e, st, sid, run in sorted(smorfs, key=lambda r: (r[0], r[1])):
        if sid in matched:
            continue
        n_new += 1
        contig_order.setdefault(c, len(contig_order))
        tag = "%s_smORF_%d" % (a.sample, n_new)
        lines.append([c, "smORF_module", "CDS", str(s), str(e), ".", st, "0",
                      "ID=%s;locus_tag=%s;product=small protein (predicted smORF);smORF=true;smorf_id=%s;"
                      "smorf_run=%s;smorf_length_aa=%d" % (tag, tag, sid, run, length.get(sid, 0))])
    lines.sort(key=lambda p: (contig_order[p[0]], int(p[3]), int(p[4])))
    with open(a.out, "w") as fo:
        fo.write("\n".join(header or ["##gff-version 3"]) + "\n")
        for p in lines:
            fo.write("\t".join(p) + "\n")
    print("%s\t%d smORFs\t%d IMP3 CDS tagged\t%d added" % (a.sample, len(smorfs), len(matched), n_new))
    assert len(matched) + n_new == len(smorfs)


if __name__ == "__main__":
    main()
