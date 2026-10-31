raw_fq_files = [
    'Preprocessing/{type}.r1.fq',
    'Preprocessing/{type}.r2.fq'
]
preprocessed_fq_files = [
    'Preprocessing/{type}.r1.preprocessed.fq.gz',
    'Preprocessing/{type}.r2.preprocessed.fq.gz',
    'Preprocessing/{type}.se.preprocessed.fq.gz',
]
trimmed_fq_files = [
    'Preprocessing/{type}.r1.trimmed.fq',
    'Preprocessing/{type}.r2.trimmed.fq',
    'Preprocessing/{type}.se.trimmed.fq'
]
rrna_filtered_fq_mt_files = [
    'Preprocessing/mt.r1.trimmed__rrna_proc.fq',
    'Preprocessing/mt.r2.trimmed__rrna_proc.fq',
    'Preprocessing/mt.se.trimmed__rrna_proc.fq',
]
bam_files = [
    'Assembly/{type}.reads.sorted.bam'
]

def reads_input_files(wildcards):
    if wildcards.type not in TYPES:
        return 'reads_input_files-no-file-here'
    if MG or MT:
        base = raw_fq_files
        if 'preprocessing' in IMP_STEPS:
            base += preprocessed_fq_files
            if PREPROCESSING_FILTERING:
                base += trimmed_fq_files
        if wildcards.type == 'mt' and 'preprocessing' in IMP_STEPS:
            base += rrna_filtered_fq_mt_files
    elif MGaln or MTaln:
        base = bam_files
    return expand(base, type=wildcards.type)

#"Analysis/mg.read_stats.tsv",

# if MG or MT:
#     rule reads_count:
#         input:
#             reads_input_files
#         output:
#             'Stats/{type}/{type}.read_stats.tsv'
#         resources:
#             runtime = VERY_SHORT_RUNTIME,
#             mem_mb = MAIN_MEMCORE_MB
#         threads:
#             1
#         # log:
#         #     "logs/analysis/read_numbers.{type}.log"
#         message:
#             "reads_count: Getting length of {wildcards.type} fastq files."
#         run:
#             for fs in {input}:
#                 for f in fs:
#                     shell('echo -e "{f}\t$(($(wc -l < {f})/4))" >> {output}')
if MGaln or MTaln:
    rule reads_count:
        input:
            reads_input_files
        output:
            'Stats/{type}/{type}.read_stats.tsv'
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            1
        conda:
            os.path.join(ENVDIR, "IMP_mapping.yaml")
        # log:
        #     "logs/analysis/read_numbers.{type}.log"
        message:
            "reads_count: Getting number of reads from {wildcards.type} bam file."
        shadow: "copy-minimal"
        shell:
             """
             echo -e "{input}\t$(samtools view -c {input})" >> {output}
             """
