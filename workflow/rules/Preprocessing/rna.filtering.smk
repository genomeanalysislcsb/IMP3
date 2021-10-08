# based on code from https://git-r3lab.uni.lu/IMP/imp3/-/tree/issue33
# see https://git-r3lab.uni.lu/IMP/imp3/-/issues/33

# https://github.com/biocore/sortmerna/
# https://github.com/biocore/sortmerna/wiki/User-manual-v4.0
rule preproc_filter_rrna:
    input:
        r1="Preprocessing/mt.r1.trimmed.fq",
        r2="Preprocessing/mt.r2.trimmed.fq",
        se="Preprocessing/mt.se.trimmed.fq",
        ref=get_rrna_refs(),
        idx=get_rrna_refs_idx(),
    output:
        filt_r1="Preprocessing/mt.r1.trimmed.rna_filtered.fq",
        filt_r2="Preprocessing/mt.r2.trimmed.rna_filtered.fq",
        filt_se="Preprocessing/mt.se.trimmed.rna_filtered.fq",
        rrna_r1="Preprocessing/mt.r1.trimmed.rna.fq",
        rrna_r2="Preprocessing/mt.r2.trimmed.rna.fq",
        rrna_se="Preprocessing/mt.se.trimmed.rna.fq",
	    otu="Preprocessing/otu_map.txt",
        # tmp files
        tmp_filt_r1=temp("Preprocessing/mt.sortmerna.filt_fwd.fq"),
        tmp_filt_r2=temp("Preprocessing/mt.sortmerna.filt_rev.fq"),
        tmp_rrna_r1=temp("Preprocessing/mt.sortmerna.rrna_fwd.fq"),
        tmp_rrna_r2=temp("Preprocessing/mt.sortmerna.rrna_rev.fq"),
    threads: 8
    resources:
        runtime="48:00:00",
        mem=get_mem_norm,
    params:
        # list of reference FASTA files
        ref=lambda wildcards, input: " ".join(["--ref %s" % ref for ref in input.ref]),
        # --idx-dir: reference index directory (w/o idx/ sub-folder)
        # --fastx: output aligned reads
        # -m: memory requirements for indexing
        params=lambda wildcards, input, threads: "--idx-dir {idx} --fastx -m {mem}".format(
            mem=get_mem_norm(wildcards, input, threads, None, "Mb"), # mem in Mb
            idx=os.path.dirname(os.path.dirname(input.idx[0]))
        ),
        # --otu-map: OTU map (input to QIIME's make_otu_table.py)
        # --paired_in: flags the PE reads as "aligned", when either of them is aligned
        # --out2: output paired reads into separate files
        params_pe="--otu_map --paired_in --out2",
        # basenames for the output files (sortmerna adds suffix/extensions automatically)
        bname_pe_filt=lambda wildcards, output: output.tmp_filt_r1.replace("_fwd.fq", ""),
        bname_pe_rrna=lambda wildcards, output: output.tmp_rrna_r1.replace("_fwd.fq", ""),
        bname_se_filt=lambda wildcards, output: os.path.splitext(output.filt_se)[0],
        bname_se_rrna=lambda wildcards, output: os.path.splitext(output.rrna_se)[0],
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log:
        "logs/preprocessing_filter_rna.log"
    message:
        "rRNA filtering: {input.r1}, {input.r2}, {input.se}"
    shell:
        """
        TMPDIR_PE=$(mktemp -d -p ${{TMPDIR}})
        TMPDIR_SE=$(mktemp -d -p ${{TMPDIR}})

        sortmerna {params.ref} {params.params} {params.params_pe} --threads {threads} \
        --reads {input.r1} --reads {input.r2} --workdir ${{TMPDIR_PE}} --aligned {params.bname_pe_rrna} --other {params.bname_pe_filt} &> {log}
        rsync -tPv {output.tmp_filt_r1} {output.filt_r1} &>> {log}
        rsync -tPv {output.tmp_filt_r2} {output.filt_r2} &>> {log}
        rsync -tPv {output.tmp_rrna_r1} {output.rrna_r1} &>> {log}
        rsync -tPv {output.tmp_rrna_r2} {output.rrna_r2} &>> {log}
        if [ ! -s {output.otu} ]; then # empty or does not exist
            touch {output.otu}
            echo "# WARNING: OTU from sortmerna file does not exists" > status/sortmerna.warning
        fi

        sortmerna {params.ref} {params.params} --threads {threads} \
        --reads {input.se} --workdir ${{TMPDIR_SE}} --aligned {params.bname_se_rrna} --other {params.bname_se_filt} &>> {log}

        cat {params.bname_pe_rrna}.log {params.bname_se_rrna}.log >> {log} && rm {params.bname_pe_rrna}.log {params.bname_se_rrna}.log
        rm -rf ${{TMPDIR_PE}} ${{TMPDIR_SE}}
        """
