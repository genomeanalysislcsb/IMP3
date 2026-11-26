#!/usr/bin/env python3
"""One GFF per annotation database from the merged annotation (annotation_CDS_RNA_hmms.gff), one line per feature
and term, for featureCounts -g <db>.

The merged annotation keeps all terms of a database in one attribute (pfam=NTP_transferase,PF00483,PF01050); counted
as is, featureCounts would see "NTP_transferase,PF00483,PF01050" as one term. Here every term gets its own line, terms
are made unique per feature, and synonyms are merged:
  pfam          Mantis reports Pfam domains by accession and by name (PF00483 and NTP_transferase); names are turned
                into accessions with the NAME / ACC pairs of Pfam-A.hmm (--pfam-hmm), versions dropped (PF00483.26)
  kegg_pathway  KEGG pathways come as koNNNNN and mapNNNNN; kept once, as koNNNNN
Without this, one domain is counted under two terms with different counts, and with featureCounts -O --fraction a
read on such a gene is split over the duplicates.
"""
import argparse, gzip, re, sys
from collections import Counter


def opener(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)


def pfam_names(path):
    """NAME -> accession (without version) from Pfam-A.hmm (or .hmm.gz / .dat)"""
    m, name = {}, None
    with opener(path) as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line.split(None, 1)[1].strip()
            elif line.startswith("ACC ") and name:
                m[name] = line.split(None, 1)[1].strip().split(".")[0]
                name = None
            elif line.startswith("#=GF ID "):                    # Pfam-A.hmm.dat
                name = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF AC ") and name:
                m[name] = line.split(None, 2)[2].strip().split(".")[0]
                name = None
    return m


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff", required=True)
    p.add_argument("--db", required=True, help="attribute of the merged GFF, e.g. kegg_ko, pfam")
    p.add_argument("--out", required=True)
    p.add_argument("--pfam-hmm", help="Pfam-A.hmm(.gz) for Pfam names -> accessions")
    a = p.parse_args()
    names = pfam_names(a.pfam_hmm) if a.db == "pfam" and a.pfam_hmm else {}
    key = re.compile(r"(?:^|;)%s=([^;]*)" % re.escape(a.db))
    n = Counter()
    with opener(a.gff) as fh, opener(a.out, "wt") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            m = key.search(f[8])
            if not m:
                continue
            terms = []
            for t in m.group(1).split(","):
                t = t.strip()
                if not t:
                    continue
                if a.db == "pfam":
                    if re.match(r"PF\d+", t):
                        t = t.split(".")[0]
                    elif t in names:
                        t, n["pfam_name_to_accession"] = names[t], n["pfam_name_to_accession"] + 1
                    else:
                        n["pfam_name_unresolved"] += 1
                elif a.db == "kegg_pathway" and t.startswith("map"):
                    t = "ko" + t[3:]
                if t not in terms:
                    terms.append(t)
                else:
                    n["duplicate_terms_dropped"] += 1
            for t in terms:
                out.write("\t".join(f[:8] + ["%s=%s" % (a.db, t)]) + "\n")
                n["lines"] += 1
            n["features"] += 1
    sys.stderr.write("split_annotation_terms %s: %s\n" % (a.db, dict(n)))


if __name__ == "__main__":
    main()
