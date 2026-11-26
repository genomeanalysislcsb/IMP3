#!/usr/bin/env python3
"""Metaproteomic validation of smORFs (smORF module of IMP3).

A smORF changes its status from 'predicted' to 'validated' when at least --min-peptides identified peptides of the
sample are unique to small proteins: the peptide occurs in the smORF and in no protein of the background (the IMP3
proteins > max_aa aa, and optionally human proteins and contaminants), as in the paper (peptides that map to both a
small protein and a protein > 100 aa were not counted). Peptides are matched as sequences, so the result does not
depend on the search database or engine (Proteome Discoverer, FragPipe, DIA-NN, ...): give the identified peptides
after FDR control.

Peptide input: a file per sample, or one table for all samples with a 'sample' column; the peptide column is
'peptide' / 'sequence' / 'Annotated Sequence' / 'Stripped.Sequence' / 'Modified.Sequence' (first match,
case-insensitive), or the first column. Modifications and flanking residues (K.PEPTIDER.G, [K].PEPTIDER.[G],
PEP(UniMod:35)TIDE, pep[+16]tide) are removed; I and L are not distinguished.
"""
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict

PEP_COLS = ["peptide", "sequence", "annotated sequence", "stripped.sequence", "modified.sequence", "peptide sequence"]
K = 6


def zopen(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def read_fasta(path):
    name, seq = None, []
    with zopen(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0] if line[1:].split() else "", []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def norm(s):
    return s.upper().replace("*", "").replace("I", "L")


def clean_peptide(p):
    p = p.strip()
    m = re.match(r"^\[?[A-Z-]*\]?\.(.+)\.\[?[A-Z-]*\]?$", p)        # K.PEPTIDER.G / [K].PEPTIDER.[G]
    if m:
        p = m.group(1)
    p = re.sub(r"^n(?=[\[(])", "", p)                                # FragPipe n[43]PEPTIDE
    p = re.sub(r"c(?=[\[(][^\])]*[\])]$)", "", p)                    # FragPipe PEPTIDEc[17]
    p = re.sub(r"\(.*?\)|\[.*?\]|\{.*?\}", "", p)                   # (UniMod:35), [+16], {..}
    p = re.sub(r"[^A-Za-z]", "", p).upper()
    return p


def read_peptides(path, sample):
    with zopen(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    if not lines:
        return set()
    delim = "\t" if "\t" in lines[0] else ","
    rows = list(csv.reader(lines, delimiter=delim))
    hdr = [h.strip().lower() for h in rows[0]]
    col = next((hdr.index(c) for c in PEP_COLS if c in hdr), None)
    scol = hdr.index("sample") if "sample" in hdr else None
    if col is None:                         # no header: first column
        col, body = 0, rows
    else:
        body = rows[1:]
    out = set()
    for r in body:
        if len(r) <= col or (scol is not None and (len(r) <= scol or r[scol].strip() != sample)):
            continue
        p = clean_peptide(r[col])
        if len(p) >= K:
            out.add(p)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--smorfs", required=True, help="smORF/<sample>.smorf.faa")
    ap.add_argument("--peptides", default="", help="identified peptides of the sample (or a table of all samples)")
    ap.add_argument("--background", nargs="*", default=[], help="FASTA of other proteins (IMP3 proteins, human, cRAP)")
    ap.add_argument("--max-aa", type=int, default=100, help="background: only proteins longer than this")
    ap.add_argument("--min-peptides", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats", required=True)
    a = ap.parse_args()

    peps = read_peptides(a.peptides, a.sample) if a.peptides else set()
    measured = bool(a.peptides) and len(peps) > 0
    by_pre = defaultdict(list)
    for p in peps:
        by_pre[norm(p)[:K]].append(norm(p))

    def scan(seq):
        """peptides (normalised) contained in seq"""
        s, hits = norm(seq), set()
        for i in range(len(s) - K + 1):
            for p in by_pre.get(s[i:i + K], ()):
                if s.startswith(p, i):
                    hits.add(p)
        return hits

    smorfs = list(read_fasta(a.smorfs))
    hits = {n: scan(s) for n, s in smorfs} if peps else {}
    found = set().union(*hits.values()) if hits else set()
    shared = set()
    if found:
        by_pre_found = defaultdict(list)
        for p in found:
            by_pre_found[p[:K]].append(p)
        by_pre = by_pre_found                 # only the peptides found in smORFs need checking
        for f in a.background:
            for n, s in read_fasta(f):
                if len(s.replace("*", "")) > a.max_aa:
                    shared |= scan(s)
    n_val = 0
    with open(a.out, "w") as fo:
        fo.write("smORF\tstatus\tmetaP_measured\tunique_peptides\tshared_peptides\tpeptides\n")
        for n, _ in smorfs:
            h = hits.get(n, set())
            uniq = sorted(h - shared)
            st = "validated" if len(uniq) >= a.min_peptides and uniq else "predicted"
            n_val += st == "validated"
            fo.write("%s\t%s\t%s\t%d\t%d\t%s\n" % (n, st, "yes" if measured else "no", len(uniq), len(h & shared),
                                                  ",".join(uniq)))
    with open(a.stats, "w") as fo:
        fo.write("sample\tmetaP_measured\tpeptides\tpeptides_in_smORFs\tpeptides_also_in_background\tsmORFs\t"
                 "validated\n")
        fo.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\n" % (a.sample, "yes" if measured else "no", len(peps), len(found),
                                                   len(shared), len(smorfs), n_val))
    sys.stderr.write("%s: %d peptides, %d smORFs validated of %d\n" % (a.sample, len(peps), n_val, len(smorfs)))


if __name__ == "__main__":
    main()
