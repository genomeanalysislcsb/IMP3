import math

rule variant_calling:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS,
        "Assembly/{type}.reads.sorted.bam"
    output:
        "Analysis/snps/{type}.variants.samtools.vcf.gz",
    resources:
        runtime = "48:00:00",
        mem = MEMCORE
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    threads:
        # Ensure at least 20gb of mem
        max(int(math.ceil(20 / (int(MEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]]))), 1)
    log:
        "logs/analysis_variant_calling.{type}.log"
    message:
        "variant_calling: Calling variants using bcftools for {wildcards.type} reads."
    shell:
        """
        #temporary directory and files
        VCF_MPU=$(mktemp --tmpdir={TMPDIR} -t "XXXXXX.mpu.vcf")
        VCF_PLT=$(mktemp --tmpdir={TMPDIR} -t "XXXXXX.plt.vcf")
        echo "$VCF_MPU" > {log}
        echo "$VCF_PLT" >> {log}

        ### run_mpileup {input[0]} {input[1]} {output[0]}
        echo "Running samtools mpileup" >> {log} 
        bcftools mpileup -A -Ou -f {input[0]} {input[1]} 2>> {log}|\
         bcftools call -vm -Ov >$VCF_MPU 2>> {log}
        bgzip -c $VCF_MPU > {output[0]} 2>> {log}
        tabix -f -p vcf {output[0]} >> {log} 2>&1
        """
