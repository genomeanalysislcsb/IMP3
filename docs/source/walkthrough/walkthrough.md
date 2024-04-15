## Run Walkthrough

This section will walk you through the test run done with the config.test.yaml, using the reads present in the test folder. It is recommended to perform this test run to check if the installation was done successfully and no steps are skipped or missed. When looking at the log, the test run will build a directed acyclic graph (DAG) of the jobs it has to run, which looks like this:

    Building DAG of jobs...
    Using shell: /usr/bin/bash
    Provided cluster nodes: 50
    Job stats:
    job                                count    min threads    max threads
    -------------------------------  -------  -------------  -------------
    ALL                                    1              1              1
    ALL_stats                              1              1              1
    ANALYSIS                               1              1              1
    ASSEMBLY                               1              1              1
    BINNING                                1              1              1
    GTDB_all                               1              1              1
    GTDBtk_summarize                       1              6              6
    PREPROCESSING                          1              1              1
    SUMMARY                                1              1              1
    TAXONOMY                               1              1              1
    allmotulinks                           1              1              1
    assembly_contig_length                 1              1              1
    bins_tar                               1              1              1
    blast_mOTU_db_genes                   10              1              1
    call_contig_depth                      1             12             12
    call_gene_depth                        2             12             12
    check_inputs                           1              1              1
    check_read_r_input                     2              1              1
    compress_binner_out                    6              8              8
    cutends                                1              1              1
    extract_mOTU_db_genes                  1              1              1
    fastqc_preprocessed                    2              3              3
    fastqc_raw                             2              2              2
    feature_count_gene                     4              8              8
    feature_count_gene_per_contig          4              8              8
    feature_count_mantis                  34              8              8
    feature_count_mantis_per_contig       34              8              8
    fetchmg                                1              6              6
    formatProteomics                       1              1              1
    gather_stats                           1              3              3
    getBins                                1              1              1
    get_present_mOTUs                      1              1              1
    index_bam                              2              1              1
    kraken_contigs                         1              5              5
    kraken_parse                           1              1              1
    kraken_reads                           2              5              5
    mapping_on_assembly                    1            128            128
    motus                                  2             18             18
    nonpareil                              2              1              1
    nonpareil_curve                        2              1              1
    parse_blast_mOTU                      10              1              1
    reads_count                            2              1              1
    reads_gzip                             2              1              1
    reads_gziped                           2              1              1
    tar_annotation                         1              8              8
    tar_intermediary                       1              1              1
    variant_calling                        2             12             12
    workflow                               1              1              1
    total                                156              1            128

