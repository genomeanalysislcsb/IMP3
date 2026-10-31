import math

rule variant_calling:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        map = "Assembly/{type}.reads.sorted.bam"
    output:
        snps = "Analysis/snps/{type}.variants.samtools.vcf.gz",
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    threads:
        # Ensure at least 20gb of mem
        max(int(math.ceil(20 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    params:
        tmp_dir=TMPDIR
    log:
        "logs/analysis/variant_calling.{type}.log"
    message:
        "variant_calling: Calling variants using bcftools for {wildcards.type} reads."
    shadow: "copy-minimal"
    shell:
        """
        # Set the temporary directory
        export TMPDIR="{params.tmp_dir}"

        # Create temporary files in the specified temporary directory
        VCF_MPU=$(mktemp "${{TMPDIR}}/XXXXXX.mpu.vcf")
        VCF_PLT=$(mktemp "${{TMPDIR}}/XXXXXX.plt.vcf")

        # Log the names of the temporary files
        echo "$VCF_MPU" > {log}
        echo "$VCF_PLT" >> {log}

        ### run_mpileup
        echo "Running samtools mpileup" >> {log} 

        # Run bcftools mpileup and call, outputting to the temporary VCF file
        bcftools mpileup -A -Ou -f {input.assembly} {input.map} 2>> {log} |\
        bcftools call -vm -Ov > "$VCF_MPU" 2>> {log}

        # Compress the VCF file and move it to the final output
        bgzip -c "$VCF_MPU" > {output.snps} 2>> {log}

        # Index the compressed VCF file
        tabix -f -p vcf {output.snps} >> {log} 2>&1
        """
