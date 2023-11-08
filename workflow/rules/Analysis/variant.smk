import math

rule variant_calling:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS,
        "Assembly/{type}.reads.sorted.bam"
    output:
        "Analysis/snps/{type}.variants.samtools.vcf.gz",
    resources:
        runtime = "24:00:00",
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    threads:
        # Ensure at least 20gb of mem
        max(int(math.ceil(20 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    params:
        tmp_dir=TMPDIR
    log:
        "logs/analysis_variant_calling.{type}.log"
    message:
        "variant_calling: Calling variants using bcftools for {wildcards.type} reads."
    shell:
        """
        #temporary directory and files
        export TMPDIR="{params.tmp_dir}"
        VCF_MPU=$(mktemp "XXXXXX.mpu.vcf")
        VCF_PLT=$(mktemp "XXXXXX.plt.vcf")
        echo "$VCF_MPU" > {log}
        echo "$VCF_PLT" >> {log}

        ### run_mpileup {input[0]} {input[1]} {output[0]}
        echo "Running samtools mpileup" >> {log} 
        bcftools mpileup -A -Ou -f {input[0]} {input[1]} 2>> {log}|\
         bcftools call -vm -Ov >$VCF_MPU 2>> {log}
        bgzip -c $VCF_MPU > {output[0]} 2>> {log}
        tabix -f -p vcf {output[0]} >> {log} 2>&1
        """
