#########################
### hmmer db counting ###
#########################

# TODO: Find a way to dynamically generate log names, since params are not accessible by log

def set_fc_runtime(wildcards, attempt):
    if attempt == 1:
        return VERY_SHORT_RUNTIME
    elif attempt == 2:
        return SHORTISH_RUNTIME
    else:
        return MEDIUM_RUNTIME

# paired total feature sum
rule feature_count_hmmer_paired:
    input:
        map = "Assembly/{type}.reads.paired.sorted.bam",
        gff = "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    output:
        temp("Analysis/annotation/{type}.hmmer.{db}.counts.paired.tsv"),
        temp("Analysis/annotation/{type}.hmmer.{db}.counts.paired.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-p -O --fraction -M --fraction",
        feature = "CDS",
        read_type = "paired",
        feature_sum = "total",
        count_type = "hmmer",
        out_dir = OUTPUTDIR,
    resources:
        runtime = set_fc_runtime,
        mem_mb = MAIN_MEMCORE_MB
    retries:
        3
    threads:
        getThreads(12)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/analysis/featurecount/hmmer/total.{type}.{db}.paired.log"
    message:
        "feature_count: Quantifying {params.db} with {wildcards.type} {params.read_type} reads ({params.feature_sum} counts)."
    shadow: "copy-minimal"
    shell:
        """
        featureCounts {params.flags} \
                      -t {params.feature} \
                      -g {params.db} \
                      -o {output[0]} \
                      -s {params.stranding} \
                      -a {input.gff} \
                      -T {threads} \
                      {input.map}> {log} 2>&1 \
                      || touch {output} \
                               {params.out_dir}/status/analysis_featureCounts.{wildcards.type}.{params.db}.{params.read_type}.impossible
        """

# single total feature sum
use rule feature_count_hmmer_paired as feature_count_hmmer_single with:
    input:
        map = "Assembly/{type}.reads.single.sorted.bam",
        gff = "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    output:
        temp("Analysis/annotation/{type}.hmmer.{db}.counts.single.tsv"),
        temp("Analysis/annotation/{type}.hmmer.{db}.counts.single.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-O --fraction -M --fraction",
        feature = "CDS",
        read_type = "single",
        feature_sum = "total",
        count_type = "hmmer",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/hmmer/total.{type}.{db}.single.log"

# paired per contig sum
use rule feature_count_hmmer_paired as feature_count_hmmer_paired_per_contig with:
    output:
        temp("Analysis/annotation/per_contig_counts/{type}.hmmer.{db}.counts.paired.tsv"),
        temp("Analysis/annotation/per_contig_counts/{type}.hmmer.{db}.counts.paired.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-p -f -O --fraction -M --fraction",
        feature = "CDS",
        read_type = "paired",
        feature_sum = "per_contig",
        count_type = "hmmer",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/hmmer/per_contig.{type}.{db}.paired.log"

# single per contig sum
use rule feature_count_hmmer_paired as feature_count_hmmer_single_per_contig with:
    input:
        map = "Assembly/{type}.reads.single.sorted.bam",
        gff = "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    output:
        temp("Analysis/annotation/per_contig_counts/{type}.hmmer.{db}.counts.single.tsv"),
        temp("Analysis/annotation/per_contig_counts/{type}.hmmer.{db}.counts.single.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-f -O --fraction -M --fraction",
        feature = "CDS",
        read_type = "single",
        feature_sum = "per_contig",
        count_type = "hmmer",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/hmmer/per_contig.{type}.{db}.single.log"

# Combine single and paired counts based on feature in first column
# (we add the 7th column of the single counts to the end of the paired counts output)
# TODO: Add summaries properly to in and outputs
rule combine_feature_counts_hmmer:
    input:
        paired = "Analysis/annotation/{type}.{db}.counts.paired.tsv",
        paired_sum = "Analysis/annotation/{type}.{db}.counts.paired.tsv.summary",
        single = "Analysis/annotation/{type}.{db}.counts.single.tsv",
        single_sum = "Analysis/annotation/{type}.{db}.counts.single.tsv.summary",
    output:
        "Analysis/annotation/{type}.hmmer.{db}.counts.tsv",
        "Analysis/annotation/{type}.hmmer.{db}.counts.tsv.summary"
    params:
        db = lambda wildcards: wildcards.db,
        feature_sum = "total",
        count_type = "hmmer",
        out_dir = OUTPUTDIR,
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log:
        "logs/analysis/featurecount/hmmer/total.{type}.{db}.combine.log"
    message:
        "feature_count: Combining {params.count_type} {params.db} {wildcards.type} paired and single read counts ({params.feature_sum} counts)."
    shadow: "copy-minimal"
    shell:
        """
        # Ensure the first columns are the same
        if [[ $(awk '{{print $1}}' {input.paired} | md5sum) != $(awk '{{print $1}}' {input.single} | md5sum) ]]; then
            echo "feature_count: The first columns of the paired and single counts are not the same." > {log}
            touch {output}
            exit 1
        fi

        # Combine the counts
        paste {input.paired} \
              <(awk -F'\t' '{{print $7}}' {input.single}) > {output[0]} 2> {log} \
              || touch {output[0]} \
                       {params.out_dir}/status/analysis_featureCounts.{wildcards.type}.{params.db}.impossible
        
        # Combine the summaries in the same way and delete the paired and single summaries
        paste {input.paired_sum} \
              <(awk -F'\t' '{{print $2}}' {input.single_sum}) > {output[1]} 2>> {log} \
              && rm {input.paired_sum} {input.single_sum} \
              || touch {output[1]}
        """

# Combine single and paired per-contig counts
use rule combine_feature_counts_hmmer as combine_feature_counts_hmmer_per_contig with:
    input:
        paired = "Analysis/annotation/per_contig_counts/{type}.{db}.counts.paired.tsv",
        paired_sum = "Analysis/annotation/per_contig_counts/{type}.{db}.counts.paired.tsv.summary",
        single = "Analysis/annotation/per_contig_counts/{type}.{db}.counts.single.tsv",
        single_sum = "Analysis/annotation/per_contig_counts/{type}.{db}.counts.single.tsv.summary"
    output:
        "Analysis/annotation/per_contig_counts/{type}.hmmer.{db}.per_contig.counts.tsv",
        "Analysis/annotation/per_contig_counts/{type}.hmmer.{db}.per_contig.counts.tsv.summary",
    params:
        db = lambda wildcards: wildcards.db,
        feature_sum = "per_contig",
        count_type = "hmmer",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/hmmer/per_contig.{type}.{db}.combine.log"

##########################
### Mantis db counting ###
##########################

# paired total feature sum
use rule feature_count_hmmer_paired as feature_count_mantis_paired with:
    input:
        map = "Assembly/{type}.reads.paired.sorted.bam",
        gff = "Analysis/annotation/intermediary/mantis/ref_gff/consensus_annotation.mantis.{db}.gff"
    output:
        temp("Analysis/annotation/{type}.mantis.{db}.counts.paired.tsv"),
        temp("Analysis/annotation/{type}.mantis.{db}.counts.paired.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-p -O --fraction -M --fraction",
        feature = "CDS",
        read_type = "paired",
        feature_sum = "total",
        count_type = "mantis",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/mantis/total.{type}.{db}.paired.log"

# single total feature sum
use rule feature_count_hmmer_paired as feature_count_mantis_single with:
    input:
        map = "Assembly/{type}.reads.single.sorted.bam",
        gff = "Analysis/annotation/intermediary/mantis/ref_gff/consensus_annotation.mantis.{db}.gff"
    output:
        temp("Analysis/annotation/{type}.mantis.{db}.counts.single.tsv"),
        temp("Analysis/annotation/{type}.mantis.{db}.counts.single.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-O --fraction -M --fraction",
        feature = "CDS",
        read_type = "single",
        feature_sum = "total",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/mantis/total.{type}.{db}.single.log"

# paired per contig sum
use rule feature_count_hmmer_paired as feature_count_mantis_paired_per_contig with:
    input:
        map = "Assembly/{type}.reads.paired.sorted.bam",
        gff = "Analysis/annotation/intermediary/mantis/ref_gff/consensus_annotation.mantis.{db}.gff"
    output:
        temp("Analysis/annotation/per_contig_counts/{type}.mantis.{db}.counts.paired.tsv"),
        temp("Analysis/annotation/per_contig_counts/{type}.mantis.{db}.counts.paired.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-p -f -O --fraction -M --fraction",
        feature = "CDS",
        read_type = "paired",
        feature_sum = "per_contig",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/mantis/per_contig.{type}.{db}.paired.log"

# single per contig sum
use rule feature_count_hmmer_paired as feature_count_mantis_single_per_contig with:
    input:
        map = "Assembly/{type}.reads.single.sorted.bam",
        gff = "Analysis/annotation/intermediary/mantis/ref_gff/consensus_annotation.mantis.{db}.gff"
    output:
        temp("Analysis/annotation/per_contig_counts/{type}.mantis.{db}.counts.single.tsv"),
        temp("Analysis/annotation/per_contig_counts/{type}.mantis.{db}.counts.single.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-f -O --fraction -M --fraction",
        feature = "CDS",
        read_type = "single",
        feature_sum = "per_contig",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/mantis/per_contig.{type}.{db}.single.log"

# Combine single and paired counts based on feature in first column
use rule combine_feature_counts_hmmer as combine_feature_counts_mantis with:
    input:
        paired = "Analysis/annotation/{type}.mantis.{db}.counts.paired.tsv",
        paired_sum = "Analysis/annotation/{type}.mantis.{db}.counts.paired.tsv.summary",
        single = "Analysis/annotation/{type}.mantis.{db}.counts.single.tsv",
        single_sum = "Analysis/annotation/{type}.mantis.{db}.counts.single.tsv.summary",
    output:
        "Analysis/annotation/{type}.mantis.{db}.counts.tsv",
        "Analysis/annotation/{type}.mantis.{db}.counts.tsv.summary"
    params:
        db = lambda wildcards: wildcards.db,
        feature_sum = "total",
        count_type = "mantis",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/mantis/total.{type}.{db}.combine.log"

# Combine single and paired per-contig counts
use rule combine_feature_counts_hmmer as combine_feature_counts_mantis_per_contig with:
    input:
        paired = "Analysis/annotation/per_contig_counts/{type}.mantis.{db}.counts.paired.tsv",
        paired_sum = "Analysis/annotation/per_contig_counts/{type}.mantis.{db}.counts.paired.tsv.summary",
        single = "Analysis/annotation/per_contig_counts/{type}.mantis.{db}.counts.single.tsv",
        single_sum = "Analysis/annotation/per_contig_counts/{type}.mantis.{db}.counts.single.tsv.summary",
    output:
        "Analysis/annotation/per_contig_counts/{type}.mantis.{db}.per_contig.counts.tsv",
        "Analysis/annotation/per_contig_counts/{type}.mantis.{db}.per_contig.counts.tsv.summary"
    params:
        db = lambda wildcards: wildcards.db,
        feature_sum = "per_contig",
        count_type = "mantis",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/mantis/per_contig.{type}.{db}.combine.log"

##########################
### Gene level counting ##
##########################

# paired total feature sum
use rule feature_count_hmmer_paired as feature_count_gene_paired with:
    input:
        map = "Assembly/{type}.reads.paired.sorted.bam",
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        temp("Analysis/annotation/{type}.{db}.{feature}_counts.paired.tsv"),
        temp("Analysis/annotation/{type}.{db}.{feature}_counts.paired.tsv.summary"),
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-p -O --fraction -M --fraction",
        feature = lambda wildcards: wildcards.feature,
        read_type = "paired",
        feature_sum = "total",
        count_type = "gene",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/gene/total.{type}.{db}.{feature}.paired.log"

# single total feature sum
use rule feature_count_hmmer_paired as feature_count_gene_single with:
    input:
        map = "Assembly/{type}.reads.single.sorted.bam",
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        temp("Analysis/annotation/{type}.{db}.{feature}_counts.single.tsv"),
        temp("Analysis/annotation/{type}.{db}.{feature}_counts.single.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-O --fraction -M --fraction",
        feature = lambda wildcards: wildcards.feature,
        read_type = "single",
        feature_sum = "total",
        count_type = "gene",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/gene/total.{type}.{db}.{feature}.single.log"

# paired per contig sum
use rule feature_count_hmmer_paired as feature_count_gene_paired_per_contig with:
    input:
        map = "Assembly/{type}.reads.paired.sorted.bam",
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        temp("Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_counts.paired.tsv"),
        temp("Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_counts.paired.tsv.summary")

    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-p -f -O --fraction -M --fraction",
        feature = lambda wildcards: wildcards.feature,
        read_type = "paired",
        feature_sum = "per_contig",
        count_type = "gene",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/gene/per_contig.{type}.{db}.{feature}.paired.log"

# single per contig sum
use rule feature_count_hmmer_paired as feature_count_gene_single_per_contig with:
    input:
        map = "Assembly/{type}.reads.single.sorted.bam",
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        temp("Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_counts.single.tsv"),
        temp("Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_counts.single.tsv.summary")
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-f -O --fraction -M --fraction",
        feature = lambda wildcards: wildcards.feature,
        read_type = "single",
        feature_sum = "per_contig",
        count_type = "gene",
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/gene/per_contig.{type}.{db}.{feature}.single.log"

# Combine single and paired counts based on feature in first column
use rule combine_feature_counts_hmmer as combine_feature_counts_gene with:
    input:
        paired = "Analysis/annotation/{type}.{db}.{feature}_counts.paired.tsv",
        paired_sum = "Analysis/annotation/{type}.{db}.{feature}_counts.paired.tsv.summary",
        single = "Analysis/annotation/{type}.{db}.{feature}_counts.single.tsv",
        single_sum = "Analysis/annotation/{type}.{db}.{feature}_counts.single.tsv.summary",
    output:
        "Analysis/annotation/{type}.{db}.{feature}_counts.tsv",
        "Analysis/annotation/{type}.{db}.{feature}_counts.tsv.summary"
    params:
        db = lambda wildcards: wildcards.db,
        feature_sum = "total",
        count_type = lambda wildcards: wildcards.feature,
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/gene/total.{type}.{db}.{feature}.combine.log"

# Combine single and paired per-contig counts
use rule combine_feature_counts_hmmer as combine_feature_counts_gene_per_contig with:
    input:
        paired = "Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_counts.paired.tsv",
        paired_sum = "Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_counts.paired.tsv.summary",
        single = "Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_counts.single.tsv",
        single_sum = "Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_counts.single.tsv.summary",
    output:
        "Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_per_contig.counts.tsv",
        "Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_per_contig.counts.tsv.summary"
    params:
        db = lambda wildcards: wildcards.db,
        feature_sum = "per_contig",
        count_type = lambda wildcards: wildcards.feature,
        out_dir = OUTPUTDIR,
    log:
        "logs/analysis/featurecount/gene/per_contig.{type}.{db}.{feature}.combine.log"
