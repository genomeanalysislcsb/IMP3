# Databases of the smORF module, downloaded into <db_path>/smorf when missing (as the paper's pipeline did)
#   Rfam covariance models (Rfam CURRENT, cmpress'ed)      ncRNA genes among the smORF candidates
#   AntiFam 8.0 HMMs                                       spurious ORFs
#   human RefSeq proteins (NCBI, deduplicated, 100 %)     metaproteomics databases (proteomics_db: human)

rule smorf_db_rfam:
    output:
        RFAM_CM
    params:
        url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
    threads: 1
    resources:
        runtime = 300,
        mem_mb = MEM
    log:
        os.path.join(OUTDIR, "logs", "db", "rfam_download.log")
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "smorf_db_rfam: downloading Rfam.cm"
    shell:
        """
        mkdir -p $(dirname {output})
        curl -sSfL {params.url} -o {output}.gz > {log} 2>&1
        curl -sSfL https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/README -o $(dirname {output})/README >> {log} 2>&1 || true
        gunzip -f {output}.gz >> {log} 2>&1
        date +%F > $(dirname {output})/DOWNLOADED
        """


rule smorf_db_rfam_press:
    input:
        RFAM_CM
    output:
        expand(RFAM_CM + "{ext}", ext=[".i1f", ".i1i", ".i1m", ".i1p"])
    threads: 1
    resources:
        runtime = 300,
        mem_mb = MEM
    log:
        os.path.join(OUTDIR, "logs", "db", "rfam_cmpress.log")
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "smorf_db_rfam_press: cmpress Rfam.cm"
    shell:
        "cmpress -F {input} > {log} 2>&1"


rule smorf_db_antifam:
    output:
        ANTIFAM_HMM
    params:
        url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/AntiFam_8.0.tar.gz"
    threads: 1
    resources:
        runtime = 60,
        mem_mb = MEM
    log:
        os.path.join(OUTDIR, "logs", "db", "antifam_download.log")
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "smorf_db_antifam: downloading AntiFam 8.0"
    shell:
        """
        d=$(dirname {output})
        mkdir -p $d
        curl -sSfL {params.url} | tar -xz -C $d > {log} 2>&1
        # AntiFam.hmm at the top of the archive (8.0) or in a sub-folder
        [[ -s {output} ]] || mv $(find $d -name AntiFam.hmm | head -n1) {output}
        echo {params.url} > $d/SOURCE
        """


rule smorf_db_human:
    output:
        HUMAN_FAA
    params:
        url = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/",
        tmp = os.path.join(TMP, "human_refseq")
    threads: T_SMALL
    resources:
        runtime = 600,
        mem_mb = MEM * T_SMALL
    log:
        os.path.join(OUTDIR, "logs", "db", "human_refseq.log")
    conda:
        os.path.join(ENVDIR, "IMP_smorf.yaml")
    message:
        "smorf_db_human: downloading and deduplicating the human RefSeq proteins"
    shell:
        """
        d=$(dirname {output})
        mkdir -p $d {params.tmp}
        files=$(curl -sSfL {params.url} | grep -o 'human\\.[0-9]*\\.protein\\.faa\\.gz' | sort -u)
        [[ -n "$files" ]] || {{ echo "no human.*.protein.faa.gz at {params.url}" > {log}; exit 1; }}
        for f in $files; do
            curl -sSfL {params.url}$f -o {params.tmp}/$f >> {log} 2>&1
        done
        cat {params.tmp}/human.*.protein.faa.gz > {params.tmp}/human_ref_final.faa.gz
        mmseqs easy-cluster {params.tmp}/human_ref_final.faa.gz {params.tmp}/human_ref_final_dedupl {params.tmp}/mm \\
            -c 1.0 --min-seq-id 1.0 --cov-mode 0 --cluster-mode 0 --threads {threads} -v 1 >> {log} 2>&1
        mv {params.tmp}/human_ref_final_dedupl_rep_seq.fasta {output}
        echo "NCBI RefSeq H_sapiens mRNA_Prot, $(date +%F): $files" > $d/SOURCE
        rm -rf {params.tmp}
        """
