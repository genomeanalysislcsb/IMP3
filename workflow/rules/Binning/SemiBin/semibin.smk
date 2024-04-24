import math

def get_min_mem_sb(wildcards, attempt):
    mem_per_core = float(config['cluster_specs']['main_cluster']['mem_per_core'])
    mem_unit = config['cluster_specs']['mem_unit']
    mem_avail_mb = math.floor(mem_per_core * to_mb_dict[mem_unit])
    return attempt * mem_avail_mb

def get_max_cores_sb(wildcards, attempt):
    mem_per_core_attempt = get_min_mem_sb(wildcards, attempt)
    max_cores = int(config['cluster_specs']['main_cluster']['cores_max'])
    total_mem = int(config['cluster_specs']['main_cluster']['mem_total'])
    attempt_supported_cores = math.floor(total_mem / mem_per_core_attempt)
    return max(min(max_cores, attempt_supported_cores), 1)

rule semibin:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        mapping = "Assembly/mg.reads.sorted.bam"
    output:
        s2b = "Binning/SemiBin/contig2bin.tsv",
        out_dir = temp(directory("Binning/SemiBin")),
        out_bin_dir = directory("Binning/SemiBin/bins")
    params:
        model = config["binning"]["SemiBin"]["model"],
        cutoff = config["binning"]["SemiBin"]["cutoff"],
        # db_dir=os.path.join(config["db_path"], "gtdb_semibin"),
        tmp_dir = os.path.join(config["tmp_dir"])
    threads:
        get_max_cores_sb
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = get_min_mem_sb
    retries:
        3
    conda:
        ENVDIR + "/IMP_semibin.yaml"
    log:
        "logs/binning/SemiBin.log"
    message:
        "Running SemiBin."
    shell:
        """
        # Run SemiBin2
        SemiBin2 \
            single_easy_bin \
                --input-fasta {input.assembly} \
                --input-bam {input.mapping} \
                --environment {params.model} \
                --threads {threads} \
                --random-seed 0 \
                --min-len {params.cutoff} \
                --compression none \
                --output {output.out_dir} >> {log} 2>&1

        # Move bin fastas and rename.
        mv {output.out_dir}/output_bins {output.out_dir}/bins
        shopt -s nullglob
        for i in {output.out_bin_dir}/*.fa; do
          bin_name=$(basename -s .fa ${{i}})
          mv ${{i}} {output.out_bin_dir}/${{bin_name}}.fasta
        done
        shopt -u nullglob

        # Write bin table
        python {SRCDIR}/write_bin_table_from_fastas.py {output.out_bin_dir} {output.s2b}
        rm -rf {output.out_dir}/output_bins
        """
