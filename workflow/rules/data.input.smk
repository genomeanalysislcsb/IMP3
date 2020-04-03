if (not 'assembly' in IMP_STEPS) and ('analysis' in IMP_STEPS or 'taxonomy' in IMP_STEPS or 'binning' in IMP_STEPS):
    if CONTIGS:
        include:
            'Assembly/common/contig.smk'
        include:
            'Assembly/common/mapping.smk'
        include:
            'Assembly/common/contig-length.smk'
        include:
            'Assembly/common/contig-depth.smk'
        if not MGaln:
            include:
                'Assembly/common/bwa.index.smk'
        if not MTaln:
            include:
                'Assembly/common/bwa.index.smk'

if (not 'analysis' in IMP_STEPS) and 'binning' in IMP_STEPS:
    if ANN:
        include:
            'Binning/DAStool/translate2DAStool.smk'
    

localrules: prepare_input_data, check_inputs, check_read_se_input, check_assembly_input, check_read_r_input

# outputs (input for workflow) to be generated (in rule)
input_files_list = []
if MG and len(MG) == 2:
    input_files_list.append("Preprocessing/mg.r1.fq")
    input_files_list.append("Preprocessing/mg.r2.fq")
elif MG and len(MG) == 3:
    input_files_list.append("Preprocessing/mg.r1.fq")
    input_files_list.append("Preprocessing/mg.r2.fq")
    input_files_list.append("Preprocessing/mg.se.fq")
elif MG and len(MG) == 1:
    input_files_list.append("Preprocessing/mg.se.fq")
if MT and len(MT) == 2:
    input_files_list.append("Preprocessing/mt.r1.fq")
    input_files_list.append("Preprocessing/mt.r2.fq")
elif MT and len(MT) == 3:
    input_files_list.append("Preprocessing/mt.r1.fq")
    input_files_list.append("Preprocessing/mt.r2.fq")
    input_files_list.append("Preprocessing/mt.se.fq")
if CONTIGS and not "assembly" in IMP_STEPS:
    input_files_list.append("Assembly/%s.assembly.merged.fa" %ASS)
    if MGaln:
        input_files_list.append("Assembly/mg.reads.sorted.bam")
    if MTaln:
        input_files_list.append("Assembly/mt.reads.sorted.bam")
if "assembly" in IMP_STEPS and LONG:
    input_files_list.append("Assembly/long_reads.fastq")
if ANN and not "analysis" in IMP_STEPS:
    input_files_list.append("Analysis/annotation/annotation_CDS_RNA_hmms.gff")

# symlinks to be generated (from , to)
symlinks_list = []
if not 'preprocessing' in IMP_STEPS and MG:
    symlinks_list.append(("mg.r1.fq", "Preprocessing/mg.r1.preprocessed.fq"))
    symlinks_list.append(("mg.r2.fq", "Preprocessing/mg.r2.preprocessed.fq"))
    symlinks_list.append(("mg.se.fq", "Preprocessing/mg.se.preprocessed.fq"))
if not 'preprocessing' in IMP_STEPS and MT:
    symlinks_list.append(("mt.r1.fq", "Preprocessing/mt.r1.preprocessed.fq"))
    symlinks_list.append(("mt.r2.fq", "Preprocessing/mt.r2.preprocessed.fq"))
    symlinks_list.append(("mt.se.fq", "Preprocessing/mt.se.preprocessed.fq"))

# symlinks to output to be generated (in rule)
symlinked_output_list = []
if not 'preprocessing' in IMP_STEPS and MG:
    symlinked_output_list.append("Preprocessing/mg.r1.preprocessed.fq")
    symlinked_output_list.append("Preprocessing/mg.r2.preprocessed.fq")
    symlinked_output_list.append("Preprocessing/mg.se.preprocessed.fq")
if not 'preprocessing' in IMP_STEPS and MT:
    symlinked_output_list.append("Preprocessing/mt.r1.preprocessed.fq")
    symlinked_output_list.append("Preprocessing/mt.r2.preprocessed.fq")
    symlinked_output_list.append("Preprocessing/mt.se.preprocessed.fq")

# empty se files to be touched
fake_output = []
if not 'preprocessing' in IMP_STEPS and MG and len(MG) == 2:
    fake_output.append("Preprocessing/mg.se.fq")
if not 'preprocessing' in IMP_STEPS and MG and len(MG) == 1:
    fake_output.append("Preprocessing/mg.r1.fq")
    fake_output.append("Preprocessing/mg.r2.fq")
if not 'preprocessing' in IMP_STEPS and MT and len(MT) == 2:
    fake_output.append("Preprocessing/mt.se.fq")


rule prepare_input_data:
    input:
        get_ori_files
    output:
        input_files_list,
        symlinked_output_list,
        fake_output
    threads: 1
    params:
        runtime = "4:00:00",
        mem = MEMCORE
    message: "Preparing input."
    run:
        prepare_input_files(input, input_files_list)
        if fake_output:
            touch_fake_output(fake_output)
        if symlinks_list:
            symlink_input_files(symlinks_list)

checked_files_list = ["status" + f.split("/")[-1] + ".checked" for f in input_files_list]

rule check_inputs:
    input:
        checked_files_list
    output:
        touch("status/inputs.done")
    
rule check_read_r_input:
    input:
        expand("Preprocessing/{wildcards.type}.{read}.fq",read=["r1","r2"])
    output:
        expand("status/{wildcards.type}.{read}.fq.checked",read=["r1","r2"])
    message: "Checking {type} read input."
    shell:
        """
        if [ "$(wc -l < {input[0]})" -eq "$(wc -l < {input[1]})" ]; then
          touch {output}
        else
          echo "ERROR: Paired read files don't have the same number of reads." && exit 1
        fi
        """

rule check_read_se_input:
    input:
        "Preprocessing/{type}.se.fq"
    output:
        "status/{type}.se.fq.checked"
    message: "Checking assembly."
    shell:
        """
        if [ -s {input} ]; then
          touch {output}
        else
          echo "ERROR: {input} is empty." && exit 1
        fi
        """
        
rule check_assembly_input:
    input:
        "Assembly/{ass}.assembly.merged.fa"
    output:
        "status/{ass}.assembly.merged.fa.checked"
    message: "Checking assembly."
    shell:
        """
        if [ -s {input} ]; then
          touch {output}
        else
          echo "ERROR: Assembly is empty." && exit 1
        fi
        """

rule check_bam_input:
    input:
        "Assembly/%s.assembly.merged.fa" %ASS,
        "Assembly/{type}.reads.sorted.bam"
    output:
        "status/{type}.reads.sorted.bam.checked"
    threads: 1
    params:
        runtime = "4:00:00",
        mem = MEMCORE
    message: "Checking {type} alignment."
    conda: ENVDIR + "/IMP_mapping.yaml"
    shell:
        """
        samtools quickcheck {input[1]} || (echo "ERROR: {input[1]} corrupt." && exit 1)
        cmp --silent <(samtools view -H {input[1]} | grep "@SQ" | cut -f 2 | sed 's/SN://' | sort) \
         <(grep ">" {input[0]} | sed 's/^>//' | sort) || (echo "ERROR: Names of assembly and {input[1]} do not match." && exit 1)
        """

rule check_ann_input:
    input:
        "Assembly/%s.assembly.merged.fa" %ASS,
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "status/annotation_CDS_RNA_hmms.gff.checked"
    threads: 1
    params:
        runtime = "4:00:00",
        mem = MEMCORE
    message: "Checking annotation."
    conda: ENVDIR + "/IMP_bioperl.yaml"
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        perl {SRCDIR}/matchContigsGFF.pl {input}
        """
