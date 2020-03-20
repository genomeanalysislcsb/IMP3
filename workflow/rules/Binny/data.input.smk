include:
    'Assembly/common/contig.rules'
include:
    'Assembly/common/mapping.rules'
include:
    '../Assembly/common/contig-length.smk'
include:
    '../Assembly/common/contig-depth.smk'
include:
    '../Assembly/common/bwa.index.smk'
    

localrules: prepare_input_data

# outputs (input for workflow) to be generated (in rule)
input_files_list = []
input_files_list.append("Assembly/%s.assembly.merged.fa" %ASS)
input_files_list.append("Assembly/mg.reads.sorted.bam")


rule prepare_input_data:
    input:
        get_ori_files
    output:
        input_files_list
    threads: 1
    params:
        runtime = "4:00:00",
        mem = "--mem-per-cpu " + MEMCORE,
        partition = "-p batch --qos qos-batch"
    message: "Preparing input."
    run:
        prepare_input_files(input, input_files_list)


