# Handle input data

# NOTE: "dummy" rule do be used for rule-inheritance
rule prepare_input_file:
    input:
        "Preprocessing/dummy.input"
    output:
        "Preprocessing/dummy.output"
    threads: 5
    resources:
        mem=get_mem_norm(5),
        runtime="01:00:00",
    message:
        "Prepare input: {input} => {output}"
    script:
        os.path.join(SRCDIR, "prepare_input_files.py")

if has_sr("mg") or has_sr("mt"):
    use rule prepare_input_file as prepare_input_file_sr with:
        input:
            lambda wildcards: config["input"]["sr"][wildcards.omic][wildcards.reads]
        output:
            "Preprocessing/{omic}.{reads}.fq"
        wildcard_constraints:
            omic="|".join(get_omics()),
            reads="|".join(get_read_types_raw("mg"))

if "assembly" in config["steps"] and has_lr("mg"):
    use rule prepare_input_file as prepare_input_file_asm with:
        input:
            config["input"]["assembly"]["contigs"]
        output:
            "Assembly/%.assembly.merged.fa" % get_omics_assembly()

    use rule prepare_input_file as prepare_input_file_asm_map with:
        input:
            lambda wildcards: config["input"]["assembly"]["mapping"][wildcards.omic]
        output:
            "Assembly/{omic}.reads.sorted.bam"

if "assembly" not in config["steps"] and has_asm():
    use rule prepare_input_file as prepare_input_file_lr_mg with:
        input:
            config["input"]["lr"]["mg"]
        output:
            "Assembly/long_reads.fastq"

if "analysis" not in config["steps"] and has_asm_gff():
    use rule prepare_input_file as prepare_input_file_asm_gff with:
        input:
            config["input"]["assembly"]["gff"]
        output:
            "Analysis/annotation/annotation_CDS_RNA_hmms.gff"

# localrules: touch_input_file
# rule touch_input_file:
#     output:
#         "Preprocessing/dummy.output"
#     message:
#         "Touch: {output}"

# if XXX:
#     use rule touch_input_file as touch_input_file_? with:
#         output:
