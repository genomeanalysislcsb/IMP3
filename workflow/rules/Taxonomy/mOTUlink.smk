COGS_LIST = config['COGS'].split()

rule allmotulinks:
    input:
        expand("Analysis/taxonomy/mOTU_links/intermediary/{COG}.bestPresentHitPhylogeny.tsv",COG=COGS_LIST)
    output:
        "Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv"
    threads:
        1
    resources:
        runtime = "02:00:00",
        mem_mb = MEMCORE_MB
    # log:
    #     "logs/analysis_mOTU_bestPresentHitPhylogeny.log"
    message:
        "allmotulinks: Extracting phylogeny for best mOTU hits."
    shell:
        """
        cat <(head -n 1 {input} | tail -n 1) <(cat {input} | grep -v queryID ) > {output}
        """

rule fetchmg:
    input:
        "Analysis/annotation/prokka.ffn",
        "Analysis/annotation/prokka.faa"
    output:
        expand("Analysis/taxonomy/mOTU_links/intermediary/marker_genes/{COG}.fna",COG=COGS_LIST)
    threads:
        getThreads(6)
    resources:
        runtime = "08:00:00",
        mem_mb = MEMCORE_MB
    params:
        outdir = "Analysis/taxonomy/mOTU_links/intermediary/marker_genes"
    log:
        "logs/analysis_fetchmg.log"
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    message:
        "fetchmg: Finding mOTU marker genes in open reading frames."
    shell:
        """
        export PATH=${{CONDA_PREFIX}}/bin:$PATH
        export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl
        {BINDIR}/fetchMG/fetchMG.pl -o {params.outdir} \
                                    -t {threads} \
                                    -d {input[0]} \
                                    -m extraction {input[1]} >> {log} 2>&1
        """

rule get_present_mOTUs:
    input:
        expand("Analysis/taxonomy/mOTUs/{type}.mOTU.counts.tsv",type=TYPES)
    output:
        "Analysis/taxonomy/mOTU_links/intermediary/all.mOTUs.present"
    threads:
        1
    resources:
        runtime = "01:00:00",
        mem_mb = MEMCORE_MB
    # log:
    #     "logs/analysis_get_present_mOTUs.log"
    message:
        "get_present_mOTUs: Extracting identity of present mOTUs."
    shell:
        """
        grep "s__" {input} | cut -f 1 | sed 's#.*s__.* \[##' | sed 's#\]##' | sort | uniq >> {output} || touch {output} status/taxonomy_mOTUs.impossible
        """

rule extract_mOTU_db_genes:
    input:
        "Analysis/taxonomy/mOTU_links/intermediary/all.mOTUs.present"
    output:
        expand("Analysis/taxonomy/mOTU_links/intermediary/{COG}.present.DBgenes.fna",COG=COGS_LIST)
    threads:
        1
    resources:
        runtime = "04:00:00",
        mem_mb = MEMCORE_MB
    params:
        outdir = "Analysis/taxonomy/mOTU_links/intermediary",
        outname = "present.DBgenes",
        coglist = ",".join(map(str, COGS_LIST)),
        # fasta = "lib/python3.9/site-packages/motus/db_mOTU/db_mOTU_DB_CEN.fasta",
        # coordinates = "lib/python3.9/site-packages/motus/db_mOTU/db_mOTU_padding_coordinates_CEN.tsv"
        fasta =  os.path.join(DBPATH,"motus_metadata","db_mOTU_DB_CEN.fasta"),
        coordinates =  os.path.join(DBPATH,"motus_metadata","db_mOTU_padding_coordinates_CEN.tsv")
    conda:
        os.path.join(ENVDIR, "IMP_motus.yaml")
    log:
        "logs/analysis_extract_mOTU_db_genes.log"
    message:
        "extract_mOTU_db_genes: Extracting genes for present mOTUs."
    shell:
        """
        export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl
        {SRCDIR}/fastaExtractCoords_mOTUs.pl -f {params.fasta} \
                                             -l {params.coordinates} \
                                             -m {input} \
                                             -c {params.coglist} \
                                             -o {params.outdir} \
                                             -n {params.outname} > {log} 2>&1 
        """

rule blast_mOTU_db_genes:
    input:
        "Analysis/taxonomy/mOTU_links/intermediary/{COG}.present.DBgenes.fna",
        "Analysis/taxonomy/mOTU_links/intermediary/marker_genes/{COG}.fna"
    output:
        "Analysis/taxonomy/mOTU_links/intermediary/{COG}.bestPresentHits.tsv",
    threads:
        1
    resources:
        runtime = "12:00:00",
        mem_mb = MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/analysis_blast_mOTU_db_genes.{COG}.log"
    message:
        "blast_mOTU_db_genes: Blasting {wildcards.COG} genes to present mOTUs."
    shell:
        """
        if [ -s {input[0]} ]; then 
          makeblastdb -in {input[0]} \
                      -dbtype nucl \
                      -parse_seqids > {log} 2>&1
          blastn -db {input[0]} \
                 -query {input[1]} \
                 -outfmt 7 \
                 -out {output} \
                 -max_target_seqs 1 >> {log} 2>&1
        else
          touch {output} status/taxonomy_mOTU_link.{wildcards.COG}.impossible
        fi
        """

rule parse_blast_mOTU:
    input:
        "Analysis/taxonomy/mOTU_links/intermediary/{COG}.bestPresentHits.tsv"
    output:
        "Analysis/taxonomy/mOTU_links/intermediary/{COG}.bestPresentHitPhylogeny.tsv"
    threads:
        1
    resources:
        runtime = "12:00:00",
        mem_mb = MEMCORE_MB
    params:
        # metaphy = "/lib/python3.9/site-packages/motus/db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv",
        # refphy = "/lib/python3.9/site-packages/motus/db_mOTU/db_mOTU_taxonomy_ref-mOTUs.tsv"
        metaphy=  os.path.join(DBPATH,"motus_metadata","db_mOTU_taxonomy_meta-mOTUs.tsv"),
        refphy =  os.path.join(DBPATH, "motus_metadata", "db_mOTU_taxonomy_ref-mOTUs.tsv")
    conda:
        ENVDIR + "/IMP_taxonomy.yaml"
    log:
        "logs/analysis_parse_blast_mOTU.{COG}.log"
    message:
        "parse_blast_mOTU: Parsing blast results of {wildcards.COG} genes."
    script:
        SRCDIR + "/mergeMOTUblast.R"

rule tar_intermediary_mOTU:
    input:
        # Make sure its done at the very end
        "status/workflow.done",
        "Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv"
    output:
        "Analysis/taxonomy/mOTU_links/intermediary.tar.gz"
    threads:
        1
    resources:
        runtime = "8:00:00",
        mem_mb = MEMCORE_MB
    params:
        intermediary = "Analysis/taxonomy/mOTU_links/intermediary/"
    log:
        "logs/mOTUlinks_tar_intermediary.log"
    message:
        "tar_intermediary_mOTU: Compressing intermediary steps of mOTU linking."
    shell:
       """
       tar cvzf {output} {params.intermediary} >> {log} 2>&1 && rm -r {params.intermediary} >> {log} 2>&1
       """



