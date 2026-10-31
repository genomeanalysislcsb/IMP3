COGS_LIST = config['COGS'].split()

rule allmotulinks:
    input:
        expand("Analysis/taxonomy/mOTU_links/intermediary/{COG}.bestPresentHitPhylogeny.tsv",COG=COGS_LIST)
    output:
        "Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv"
    threads:
        1
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    # log:
    #     "logs/analysis/mOTU_bestPresentHitPhylogeny.log"
    message:
        "allmotulinks: Extracting phylogeny for best mOTU hits."
    shadow: "copy-minimal"
    shell:
        """
        cat <(head -n 1 {input} | tail -n 1) <(cat {input} | grep -v queryID ) > {output}
        """

rule fetchmg:
    input:
        "Analysis/annotation/bakta/annotation.ffn",
        "Analysis/annotation/bakta/annotation.faa"
    output:
        expand("Analysis/taxonomy/mOTU_links/intermediary/marker_genes/{COG}.fna",COG=COGS_LIST)
    threads:
        getThreads(6)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        outdir = "Analysis/taxonomy/mOTU_links/intermediary/marker_genes"
    log:
        "logs/analysis/fetchmg.log"
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    message:
        "fetchmg: Finding mOTU marker genes in open reading frames."
    shadow: "copy-minimal"
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
        expand("Analysis/taxonomy/mOTUs/{type}.mOTU.{out_type}.tsv", type=TYPES, out_type=["counts", "rel_abun"])
    output:
        "Analysis/taxonomy/mOTU_links/intermediary/all.mOTUs.present"
    params:
        impossible_check = f'{OUTPUTDIR}/status/taxonomy_mOTUs.impossible'
    threads:
        1
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    # log:
    #     "logs/analysis/get_present_mOTUs.log"
    message:
        "get_present_mOTUs: Extracting identity of present mOTUs."
    shadow: "copy-minimal"
    shell:
        """
        grep "s__" {input} | cut -f 1 | sed 's#.*s__.* \[##' | sed 's#\]##' | sort | uniq >> {output} || touch {output} {params.impossible_check}
        """

rule extract_mOTU_db_genes:
    input:
        "Analysis/taxonomy/mOTU_links/intermediary/all.mOTUs.present"
    output:
        expand("Analysis/taxonomy/mOTU_links/intermediary/{COG}.present.DBgenes.fna",COG=COGS_LIST)
    threads:
        1
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
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
        "logs/analysis/taxonomy/extract_mOTU_db_genes.log"
    message:
        "extract_mOTU_db_genes: Extracting genes for present mOTUs."
    shadow: "copy-minimal"
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
    params:
        out_dir = OUTPUTDIR
    threads:
        1
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/analysis/taxonomy/post_process_motus/blast_mOTU_db_genes.{COG}.log"
    message:
        "blast_mOTU_db_genes: Blasting {wildcards.COG} genes to present mOTUs."
    shadow: "copy-minimal"
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
          touch {output} {params.out_dir}/status/taxonomy_mOTU_link.{wildcards.COG}.impossible
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
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        # metaphy = "/lib/python3.9/site-packages/motus/db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv",
        # refphy = "/lib/python3.9/site-packages/motus/db_mOTU/db_mOTU_taxonomy_ref-mOTUs.tsv"
        metaphy=  os.path.join(DBPATH,"motus_metadata","db_mOTU_taxonomy_meta-mOTUs.tsv"),
        refphy =  os.path.join(DBPATH, "motus_metadata", "db_mOTU_taxonomy_ref-mOTUs.tsv")
    conda:
        ENVDIR + "/IMP_taxonomy.yaml"
    log:
        "logs/analysis/taxonomy/post_process_motus/parse_blast_mOTU.{COG}.log"
    message:
        "parse_blast_mOTU: Parsing blast results of {wildcards.COG} genes."
    script:
        SRCDIR + "/mergeMOTUblast.R"

rule tar_intermediary_mOTU:
    input:
        # Make sure its done at the very end
        check = "status/workflow.done",
        file = "Analysis/taxonomy/mOTU_links/bestPresentHitPhylogeny.tsv"
    output:
        "Analysis/taxonomy/mOTU_links/intermediary.tar.gz"
    threads:
        1
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        intermediary = "Analysis/taxonomy/mOTU_links/intermediary/"
    log:
        "logs/misc/mOTUlinks_tar_intermediary.log"
    message:
        "tar_intermediary_mOTU: Compressing intermediary steps of mOTU linking."
    shadow: "copy-minimal"
    shell:
       """
       tar cvzf {output} {params.intermediary} >> {log} 2>&1 && rm -r {params.intermediary} >> {log} 2>&1
       """



