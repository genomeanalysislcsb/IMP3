rule cp_seqs:
    input:
        lambda wildcards: expand("{raw_directory}/Analysis/annotation/prokka.f{{ext}}",
                raw_directory=samples.path[wildcards.sample])
    output:
        temporary("Catalogue/{sample}.prokka.all.f{ext}")
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
#    log: "logs/catalogue_cp.{sample}.{ext}.log"
    message: 'catalogue_cp_seqs: Making temporary copy of sequences.'
    threads: 1
    shell:
        """
        cp {input} {output}
        """

rule get_complete_proteins:
    input:
        lambda wildcards: expand("{raw_directory}/Analysis/annotation/annotation_CDS_RNA_hmms.gff",
                raw_directory=samples.path[wildcards.sample]),
        "Catalogue/{sample}.prokka.all.f{ext}"
    output:
        temporary("Catalogue/{sample}.complete_proteins.f{ext}"),
        temporary("Catalogue/{sample}.prokka.all.f{ext}.index")
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_bioperl.yaml"
    log: "logs/catalogue_get_complete_proteins.{sample}.f{ext}.log"
    message: 'catalogue_get_complete_proteins: Extracting complete CDS from {wildcards.sample} prokka sequences.'
    threads: 1
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2  
        {SRCDIR}/fastaExtract_complete.pl {input[1]} {input[0]} >> {output[0]}
        """

rule cat:
    input:
        expand("Catalogue/{sample}.complete_proteins.f{{ext}}",
                sample=samples.loc[samples.genes == 1,"sample"])
    output:
        "Catalogue/multi_prokka.f{ext}"
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
#    log: "logs/catalogue_cat.log"
    message: 'catalogue_cat: concatenating complete coding regions.'
    threads: 1
    shell:
        """
        cat {input} >> {output}
        """

rule cdhit:
    input:
        "Catalogue/multi_prokka.ffn"
    output:
        "Catalogue/catalogue.fasta"
    threads: getThreads(BIGCORENO)
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    params:
        tmp = TMPDIR,
        cat = "Catalogue/catalogue",
        limmem = BIGMEMTOTAL * 1000 
    log: "logs/catalogue_cdhit.log"
    message: 'catalogue_cdhit: clustering coding regions.'
    conda: ENVDIR + "/IMP_multi.yaml"
    shell:
        """
        cd-hit-est -c 0.95 -T {threads} -M {params.limmem} -G 0 -aS 0.9 -g 1 -r 0 -d 0 -n 10 -i {input} -o {params.cat} >> {log} 2>&1
        mv {params.cat} {output}
        
        """

rule mmSeqs:
    input:
        "Catalogue/multi_prokka.faa"
    output:
        "Catalogue/aa_catalogue_cluster.tsv",
        "Catalogue/aa_catalogue.faa"
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    params:
        tmp = TMPDIR,
        cat = "Catalogue/aa_catalogue"
    log: "logs/catalogue_mmSeqs.log"
    message: 'catalogue_mmSeqs: clustering amino acid sequences.'
    conda: ENVDIR + "/IMP_multi.yaml"
    shell:
        """
        mmseqs easy-cluster {input} {params.cat} {params.tmp} --min-seq-id 0.95 -c 0.9 --cov-mode 0 >> {log} 2>&1
        mv {params.cat}_rep_seq.fasta {output[1]}
        """
       
rule extract_gff:
    input:
        "Catalogue/{sample}.rep.names",
        lambda wildcards: expand("{raw_directory}/Analysis/annotation/annotation_CDS_RNA_hmms.gff",
                raw_directory=samples.loc[wildcards.sample ,"path"])
    output:
        temporary("Catalogue/catalogue.{sample}.gff")
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        tmp = TMPDIR
    log: "logs/catalogue_extract_gff.{sample}.log"
    message: 'collate_extract_gff: getting annotations for representative genes from {wildcards.sample} gff.'
    threads: 1
    script:
        SRCDIR + "/catalogue_extract_gff.py"
       
rule cat_gff:
    input:
        expand("Catalogue/catalogue.{sample}.gff",
                sample=samples.loc[samples.genes == 1 ,"sample"])
    output:
        temporary("Catalogue/catalogue.gff")
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        tmp = TMPDIR
    log: "logs/catalogue_cat_gff.log"
    message: 'collate_cat_gff: concatenating annotations for representative genes.'
    threads: 1
    shell:
        """
        cat {input} >> {output}
        """

rule extract_representative:
    input:
        lambda wildcards: expand("{raw_directory}/Analysis/annotation/prokka.faa",
                raw_directory=samples.path[wildcards.sample]),
        "Catalogue/catalogue.fasta"
    output:
        temporary('Catalogue/{sample}.prokka.rep.faa'),
        temporary('Catalogue/{sample}.prokka.rep.faa.index'),
        temporary('Catalogue/{sample}.rep.names'),
        temporary('Catalogue/{sample}.rep.faa')
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        prok_p = lambda wildcards: samples.prokka_prefix[wildcards.sample] + "_"
    threads: 1
    conda: ENVDIR + "/IMP_bioperl.yaml"
    log: "logs/catalogue_extract_representative.{sample}.log"
    message: "extract_representative: Extracting representative amino acid sequences from {wildcards.sample}."
    shell:
        """
        cp {input[0]} {output[0]}
        sed '/^>/s/>//' {input[1]} | sed 's/ .+//' | grep {params.prok_p} >> {output[2]}
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SRCDIR}/fastaExtract.pl {output[0]} {output[2]} >> {output[3]}
        """

rule concat_reprentatitves:
    input:
        expand("Catalogue/{sample}.prokka.rep.faa",
                sample=samples.loc[samples.prokka_prefix == samples.prokka_prefix,"sample"])
    output:
        'Catalogue/catalogue.faa'
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    threads: 1
    log: "logs/catalogue_concat_representatives.log"
    message: "concat_representatives: Concatenating all representative amino acid sequences."
    shell:
        """
        cat {input} >> {output}
        """


