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
    elif 'preprocessing' in IMP_STEPS:
        # issue #28: required for new FASTAs to be used for filtering
        include:
            'Assembly/common/bwa.index.smk'
    
localrules: check_inputs, check_read_se_input, check_assembly_input, check_longread_input

# outputs (input for workflow) to be generated (in rule)
input_files_list = []
if MG and len(MG) == 2:
    input_files_list.append("Preprocessing/mg.r1.fq.gz")
    input_files_list.append("Preprocessing/mg.r2.fq.gz")
elif MG and len(MG) == 3:
    input_files_list.append("Preprocessing/mg.r1.fq.gz")
    input_files_list.append("Preprocessing/mg.r2.fq.gz")
    input_files_list.append("Preprocessing/mg.se.fq.gz")
elif MG and len(MG) == 1:
    input_files_list.append("Preprocessing/mg.se.fq.gz")
if MT and len(MT) == 2:
    input_files_list.append("Preprocessing/mt.r1.fq.gz")
    input_files_list.append("Preprocessing/mt.r2.fq.gz")
elif MT and len(MT) == 3:
    input_files_list.append("Preprocessing/mt.r1.fq.gz")
    input_files_list.append("Preprocessing/mt.r2.fq.gz")
    input_files_list.append("Preprocessing/mt.se.fq.gz")
elif MT and len(MT) == 1:
    input_files_list.append("Preprocessing/mt.se.fq.gz")
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
    symlinks_list.append(("mg.r1.fq.gz", "Preprocessing/mg.r1.preprocessed.fq.gz"))
    symlinks_list.append(("mg.r2.fq.gz", "Preprocessing/mg.r2.preprocessed.fq.gz"))
    symlinks_list.append(("mg.se.fq.gz", "Preprocessing/mg.se.preprocessed.fq.gz"))
if not 'preprocessing' in IMP_STEPS and MT:
    symlinks_list.append(("mt.r1.fq.gz", "Preprocessing/mt.r1.preprocessed.fq.gz"))
    symlinks_list.append(("mt.r2.fq.gz", "Preprocessing/mt.r2.preprocessed.fq.gz"))
    symlinks_list.append(("mt.se.fq.gz", "Preprocessing/mt.se.preprocessed.fq.gz"))

# symlinks to output to be generated (in rule)
symlinked_output_list = []
if not 'preprocessing' in IMP_STEPS and MG:
    symlinked_output_list.append("Preprocessing/mg.r1.preprocessed.fq.gz")
    symlinked_output_list.append("Preprocessing/mg.r2.preprocessed.fq.gz")
    symlinked_output_list.append("Preprocessing/mg.se.preprocessed.fq.gz")
if not 'preprocessing' in IMP_STEPS and MT:
    symlinked_output_list.append("Preprocessing/mt.r1.preprocessed.fq.gz")
    symlinked_output_list.append("Preprocessing/mt.r2.preprocessed.fq.gz")
    symlinked_output_list.append("Preprocessing/mt.se.preprocessed.fq.gz")

# empty se files to be touched
fake_output = []
if not 'preprocessing' in IMP_STEPS and MG and len(MG) == 2:
    fake_output.append("Preprocessing/mg.se.fq.gz")
if MG and len(MG) == 1:
    fake_output.append("Preprocessing/mg.r1.fq.gz")
    fake_output.append("Preprocessing/mg.r2.fq.gz")
if not 'preprocessing' in IMP_STEPS and MT and len(MT) == 2:
    fake_output.append("Preprocessing/mt.se.fq.gz")
if MT and len(MT) == 1:
    fake_output.append("Preprocessing/mt.r1.fq.gz")
    fake_output.append("Preprocessing/mt.r2.fq.gz")

localrules: prepare_input_data_touch

rule prepare_input_data_touch:
    input:
        get_ori_files
    output:
        fake_output
    message:
        "Prepare input: touch"
    run:
        if fake_output:
            touch_fake_output(fake_output)

localrules: prepare_input_data_symlink

rule prepare_input_data_symlink:
    input:
       expand([
        'Preprocessing/{type}.r1.fq.gz',
        'Preprocessing/{type}.r2.fq.gz',
        'Preprocessing/{type}.se.fq.gz'], type=TYPES) 
    output:
        symlinked_output_list
    message:
        "Prepare input: symlinks"
    run:
        if symlinks_list:
            symlink_input_files(symlinks_list)

rule prepare_input_data_copy:
    input:
        get_ori_files
    output:
        temp(input_files_list),
    threads:
        getThreads(5)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "Prepare input: decompress/copy"
    script:
        os.path.join(SRCDIR, "prepare_input_files.py")

checked_files_list = ["status/" + f.split("/")[-1] + ".checked" for f in input_files_list]

rule check_inputs:
    input:
        checked_files_list
    output:
        touch("status/inputs.done")
    priority: 11
    
rule check_read_r_input:
    input:
        r1="Preprocessing/{type}.r1.fq.gz",
        r2="Preprocessing/{type}.r2.fq.gz"
    output:
        expand("status/{{type}}.{read}.fq.gz.checked",read=["r1","r2"]) 
    priority: 10
    threads: 8
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log:
        'logs/check_inputs_{type}_pe.log'
    message:
        "Checking paired read counts for {wildcards.type} using seqkit."
    shell:
        """
        # Use seqkit stats to get read counts efficiently
        # -j {threads} uses multiple threads
        # -T produces tab-separated output without headers
        # cut -f 5 extracts the 'num_seqs' column (usually the 5th field)

        seqkit_exit_code=0
        stats=$(seqkit stats -j {threads} -T {input.r1} {input.r2} 2>> {log}) || seqkit_exit_code=$?

        if [ $seqkit_exit_code -ne 0 ]; then
            echo "seqkit stats command failed with exit code $seqkit_exit_code for {wildcards.type}. Check log {log}" >&2
            exit $seqkit_exit_code
        fi

        # Extract the 4th field (num_seqs) from the 2nd line (R1 stats) and 3rd line (R2 stats) of the captured output.
        # Assumes seqkit succeeded and produced a header + two data lines.
        count_r1=$(awk 'NR==2 {{print $4}}' <<< "$stats")
        count_r2=$(awk 'NR==3 {{print $4}}' <<< "$stats")

        # Check if counts were obtained (non-empty), are greater than zero, and are equal
        if [[ -n "$count_r1" && "$count_r1" -gt 0 && "$count_r1" -eq "$count_r2" ]]; then
            # If counts match and are non-zero, touch the output file to signal success
            touch {output} 2>> {log}
        else
            # If counts differ, are zero, or couldn't be determined, report error and exit
            # Redirecting echo to stderr (>&2) is good practice for error messages
            echo "R1 ({input.r1}): $count_r1 reads"  &>> {log}
            echo "R2 ({input.r2}): $count_r2 reads"  &>> {log}
            exit 1 # Ensure Snakemake knows the rule failed
        fi
        """


rule check_read_se_input:
    input:
        "Preprocessing/{type}.se.fq.gz"
    output:
        "status/{type}.se.fq.gz.checked"
    priority: 10
    threads: 4
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log:
        "logs/check_inputs_{type}_se.log"
    message:
        "Checking single-end reads in {input} using seqkit (capture/parse)."
    shell:
        """
        # Use seqkit stats to count reads in the gzipped FASTQ
        # Capture stdout, check exit code, then parse. Redirect stderr to log.

        seqkit_exit_code=0
        stats=$(seqkit stats -j {threads} -T {input} 2>> {log}) || seqkit_exit_code=$?

        if [ $seqkit_exit_code -ne 0 ]; then
            echo "seqkit stats command failed with exit code $seqkit_exit_code for {wildcards.type} SE input {input}. Check log {log}" >&2
            exit $seqkit_exit_code
        fi

        # Extract the 4th field (num_seqs) from the 2nd line (the first data line after header)
        # Using awk on the variable via a here-string <<<"$stats"
        count_se=$(awk 'NR==2 {{print $4}}' <<< "$stats")

        # Check if count was obtained (non-empty) and is greater than zero
        if [[ -n "$count_se" && "$count_se" -gt 0 ]]; then
            # Success: file is readable and contains reads
            # Create the empty output file, redirect potential touch errors to log
            touch {output} 2>> {log}
        else
            # Failure: file contains zero reads, is unreadable, or parsing failed
            # Log detailed info before exiting. Using &>> appends both stdout and stderr of echo.
            echo "Read count check failed for SE input {input}." &>> {log}
            echo "Reason: Count zero, negative, or extraction failed." &>> {log}
            echo "--- Seqkit Stats Output ---" &>> {log}
            echo "$stats" &>> {log}
            echo "--- Extracted Count ---" &>> {log}
            echo "SE ({input}): '$count_se' reads" &>> {log}
            # Also print a concise error to stderr for Snakemake console output
            echo "Read count check failed for SE input {input}. Count=\'$count_se\'. See log {log}" >&2
            exit 1 # Ensure Snakemake knows the rule failed
        fi
        """
        
rule check_assembly_input:
    input:
        "Assembly/{ass}.assembly.merged.fa"
    output:
        "status/{ass}.assembly.merged.fa.checked"
    priority: 10
    threads: 1
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log:
        "logs/check_inputs_{ass}.log"
    message:
        "Checking assembly sequences in {input} using seqkit (capture/parse)."
    shell:
        r"""
        # Use seqkit stats to count sequences in the FASTA file.
        # NOTE: -T produces tab-separated output *including* a header line.
        # Capture stdout, check exit code, then parse. Redirect stderr to log.

        seqkit_exit_code=0
        stats=$(seqkit stats -j {threads} -T {input} 2>> {log}) || seqkit_exit_code=$?

        if [ $seqkit_exit_code -ne 0 ]; then
            echo "seqkit stats command failed with exit code $seqkit_exit_code for assembly {input}. Check log {log}" >&2
            exit $seqkit_exit_code
        fi

        # Extract the 4th field (num_seqs) from the 2nd line (the first data line after header)
        # Using awk on the variable via a here-string <<<"$stats"
        count_ass=$(awk 'NR==2 {{print $4}}' <<< "$stats")

        # Check if count was obtained (non-empty) and is greater than zero
        if [[ -n "$count_ass" && "$count_ass" -gt 0 ]]; then
            # Success: file is readable and contains sequences
            # Create the empty output file, redirect potential touch errors to log
            touch {output} 2>> {log}
        else
            # Failure: file contains zero sequences, is unreadable, or parsing failed
            # Log detailed info before exiting. Using &>> appends both stdout and stderr of echo.
            echo "Sequence count check failed for assembly {input}." &>> {log}
            echo "Reason: Count zero, negative, or extraction failed." &>> {log}
            echo "--- Seqkit Stats Output ---" &>> {log}
            echo "$stats" &>> {log}
            echo "--- Extracted Count ---" &>> {log}
            echo "Assembly ({input}): '$count_ass' sequences" &>> {log}
            # Also print a concise error to stderr for Snakemake console output
            echo "Sequence count check failed for assembly {input}. Count='${count_ass}'. See log {log}" >&2
            exit 1 # Ensure Snakemake knows the rule failed
        fi
        """

rule check_longread_input:
    input:
        "Assembly/long_reads.fastq"
    output:
        "status/long_reads.fastq.checked"
    priority: 10
    message: "Checking long reads."
    shell:
        """
        if [ -s {input} ]; then
          touch {output}
        else
          echo "ERROR: Long reads file is empty." && exit 1
        fi
        """

rule check_bam_input:
    input:
        "Assembly/%s.assembly.merged.fa" %ASS,
        "Assembly/{type}.reads.sorted.bam"
    output:
        "status/{type}.reads.sorted.bam.checked"
    priority: 10
    threads: 1
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    message: "Checking {wildcards.type} alignment."
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    shadow: "copy-minimal"
    shell:
        """
        samtools quickcheck {input[1]} || (echo "ERROR: {input[1]} corrupt." && exit 1)
        (cmp --silent <(samtools view -H {input[1]} | grep "@SQ" | cut -f 2 | sed 's/SN://' | sort) \
         <(grep ">" {input[0]} | sed 's/^>//' | sort) && touch {output} )|| (echo "ERROR: Names of assembly and {input[1]} do not match." && exit 1)
        """

rule check_ann_input:
    input:
        "Assembly/%s.assembly.merged.fa" %ASS,
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "status/annotation_CDS_RNA_hmms.gff.checked"
    priority: 10
    threads: 1
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    message: "Checking annotation."
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")   # perl-bioperl (Bio::DB::Fasta)
    shadow: "copy-minimal"
    shell:
        """
        if [ -s {input[1]} ]; then
          export PATH=${{CONDA_PREFIX}}/bin:$PATH
          export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl
          (perl {SRCDIR}/matchContigsGFF.pl {input} && touch {output} ) || (echo "ERROR: Assembly and {input[1]} do not match." && exit 1)
        else
          echo "ERROR: Annotation file is empty." && exit 1
        fi
        """
