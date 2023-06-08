if MG or MT:
    rule reads_count_sum:
        input:
            reads_input_files
        output:
            'Stats/{type}/{type}.read_counts.txt'
        resources:
            runtime = "04:00:00",
            mem_mb = MAIN_MEMCORE_MB
        threads:
            1
        # log:
        #     "logs/analysis_read_numbers.{type}.log"
        message:
            "reads_count: Getting length of {wildcards.type} fastq files."
        run:
            for fs in {input}:
                for f in fs:
                    shell('echo -e "{f}\t$(($(wc -l < {f})/4))" >> {output}')
elif MGaln or MTaln:
    rule reads_count_sum:
        input:
            reads_input_files
        output:
            'Stats/{type}/{type}.read_counts.txt'
        resources:
            runtime = "04:00:00",
            mem_mb = MAIN_MEMCORE_MB
        threads:
            1
        conda:
            os.path.join(ENVDIR, "IMP_mapping.yaml")
        # log:
        #     "logs/analysis_read_numbers.{type}.log"
        message:
            "reads_count: Getting number of reads from {wildcards.type} bam file."
        shell:
             """
             echo -e "{input}\t$(samtools view -c {input})" >> {output}
             """
