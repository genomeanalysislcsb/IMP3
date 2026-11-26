#!/usr/bin/env python3
"""Headers of the metaproteomics search databases (smORF module of IMP3).

Same header structure as the per-sample databases of the smORF methods paper (MP_DB: copy_bakta*.py,
modify_headers_uniqueID.py, add_smorf_tag.py, get_smorfs_with_new_headers_diff_IDs.py, cRAP / HumanRef
obtain_headers.py by Luis F. Delgado and Julia Ortis Sunyer, GPL-3.0), in one script:

  imp    Bakta proteins of an IMP3 run      >imp|ID=<prefix>_imp_<sample>_<n>|gene=<header>|contig=<contig>|sample=<sample>|smORF=<True|False>
  smorf  smORFs of the prediction            >smorf|ID=<prefix>_smorf_<sample>_<n>|gene=<gene>|contig=<contig>|sample=<sample>|smORF=True
  human  human reference proteome           >human|ID=<prefix>_human_<n>|gene=<header>
  crap   common contaminants (CCP cRAP)     >crap|ID=<prefix>_cRAP_<nnn>|gene=<accession> <description>

Differences to the paper's scripts: the protein number is counted per sample and the ID contains the sample (the
paper numbered all samples together), so a sample can be added later without renumbering the others; white space
in headers is replaced by '_' (MMseqs2 cuts names at the first white space, and the cluster IDs are parsed from the
full header).
"""
import argparse
import gzip
import re
import sys


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
                name, seq = line[1:].strip(), []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def clean(s):
    """no white space and no '|' inside a field"""
    return re.sub(r"\s+", "_", s.strip()).replace("|", "/")


def bakta_contigs(tsv):
    """locus tag -> contig from Bakta's annotation.tsv (Sequence Id, Type, Start, Stop, Strand, Locus Tag, ...)"""
    out = {}
    with zopen(tsv) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                p = line.split()
            if len(p) >= 6:
                out[p[5]] = p[0]
    return out


def write(records, out, width):
    with open(out, "w") as fh:
        for header, seq in records:
            fh.write(">%s\n%s\n" % (header, seq.replace("*", "")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["imp", "smorf", "human", "crap"])
    ap.add_argument("--fasta", required=True, help="input proteins (.faa, gzip or not)")
    ap.add_argument("--tsv", help="imp: Bakta annotation.tsv (locus tag -> contig)")
    ap.add_argument("--sample", default="", help="imp, smorf: sample name")
    ap.add_argument("--prefix", default="imp3", help="first part of the protein IDs (paper: expobiome)")
    ap.add_argument("--max-len", type=int, default=100, help="imp: smORF=True up to this length (aa)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    recs = list(read_fasta(a.fasta))
    width = len(str(max(len(recs), 1)))
    out = []
    if a.kind == "imp":
        if not a.tsv:
            sys.exit("imp needs --tsv")
        contigs = bakta_contigs(a.tsv)
        for i, (h, s) in enumerate(recs, 1):
            gene = h.split()[0]
            tag = "True" if len(s.replace("*", "")) <= a.max_len else "False"
            out.append(("imp|ID=%s_imp_%s_%s|gene=%s|contig=%s|sample=%s|smORF=%s"
                        % (a.prefix, a.sample, str(i).zfill(width), clean(h), contigs.get(gene, "NA"), a.sample, tag), s))
    elif a.kind == "smorf":
        for i, (h, s) in enumerate(recs, 1):
            gene = h.split("::", 1)[-1].split("#")[0].strip()          # default_<contig>_<n> / small_<contig>_<n>
            contig = gene.split("_", 1)[-1].rsplit("_", 1)[0]
            out.append(("smorf|ID=%s_smorf_%s_%s|gene=%s|contig=%s|sample=%s|smORF=True"
                        % (a.prefix, a.sample, str(i).zfill(width), clean(gene), contig, a.sample), s))
    elif a.kind == "human":
        for i, (h, s) in enumerate(recs, 1):
            out.append(("human|ID=%s_human_%s|gene=%s" % (a.prefix, str(i).zfill(width), clean(h)), s))
    else:   # crap: CCP cRAP (>sp|cRAP001|P00330|ADH1_YEAST ...) or already converted (>crap|ID=..._cRAP_001|gene=...)
        for i, (h, s) in enumerate(recs, 1):
            m = re.match(r"crap\|ID=\S*?_cRAP_(\d+)\|gene=(.*)$", h)
            if m:
                num, rest = m.group(1), m.group(2)
            else:
                p = h.split("|")
                if len(p) >= 4:
                    mm = re.match(r"cRAP(\d+)", p[1])
                    num = mm.group(1).zfill(3) if mm else str(i).zfill(3)
                    rest = "%s %s" % (p[2], "|".join(p[3:]))
                else:
                    num, rest = str(i).zfill(3), h
            out.append(("crap|ID=%s_cRAP_%s|gene=%s" % (a.prefix, num, clean(rest)), s))
    write(out, a.out, width)
    sys.stderr.write("%s: %d proteins -> %s\n" % (a.kind, len(out), a.out))


if __name__ == "__main__":
    main()
