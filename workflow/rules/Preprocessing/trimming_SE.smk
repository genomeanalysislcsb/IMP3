import os

# TODO: CHeck what should take predecence: adapter_fasta or adapter_sequence
def get_fastp_adapter_params_se(wildcards):
    adapter_arg = ""
    # Check if fasta is defined in config and use it if it is
    adapter_fasta = config["fastp"]["adapter_trimming"][wildcards.type]["adapter_sequence_fasta"]
    if adapter_fasta:
        return f"--adapter_fasta {adapter_fasta}"

    # Get adapter sequence from config
    r1_adapt_str = config["fastp"]["adapter_trimming"][wildcards.type]["adapter_sequence_r1"]

    # If adapter_sequence_r1 is defined in config, use it
    if r1_adapt_str:
        adapter_arg = f"--adapter_sequence {r1_adapt_str}"

    return adapter_arg

def get_fastp_params(wildcards):
    # Adapter trimming
    params = [get_fastp_adapter_params_se(wildcards)]

    # Quality filtering
    qf = config["fastp"]["quality_filtering"]
    params.extend([
        f"--qualified_quality_phred {qf['qualified_quality_phred']}",
        f"--unqualified_percent_limit {qf['unqualified_percent_limit']}",
        f"--n_base_limit {qf['n_base_limit']}",
        f"--average_qual {qf['average_qual']}",
        "--low_complexity_filter" if qf.get("low_complexity_filter", False) else "",
        f"--complexity_threshold {qf['complexity_threshold']}"
    ])

    # PolyG and PolyX tail trimming
    if config["fastp"]["polyG_tail_trimming"]["force_trim_poly_g"]:
        params.append(f"--trim_poly_g --poly_g_min_len {config['fastp']['polyG_tail_trimming']['poly_g_min_len']}")
    elif config["fastp"]["polyG_tail_trimming"]["force_disable_trim_poly_g"]:
        params.append("--disable_trim_poly_g")

    if config["fastp"]["polyX_tail_trimming"]["trim_poly_x"]:
        params.append(f"--trim_poly_x --poly_x_min_len {config['fastp']['polyX_tail_trimming']['poly_x_min_len']}")

    # Global trimming
    trim_front1 = config["fastp"]["global_trimming"]["trim_front1"]
    trim_front2 = config["fastp"]["global_trimming"]["trim_front2"]
    trim_tail = config["fastp"]["global_trimming"]["trim_tail1"]
    trim_tail2 = config["fastp"]["global_trimming"]["trim_tail2"]

    if trim_front2:
        params.extend([
            f"--trim_front1 {trim_front1}",
            f"--trim_front2 {trim_front2}"
        ])
    else:
        params.append(f"--trim_front1 {trim_front1}")
    
    if trim_tail2:
        params.extend([
            f"--trim_tail1 {trim_tail}",
            f"--trim_tail2 {trim_tail2}"
        ])
    else:
        params.append(f"--trim_tail1 {trim_tail}")

    # Quality cutting
    qc = config["fastp"]["quality_cutting"]
    if qc["cut_front"]:
        params.extend([
            f"--cut_front --cut_front_window_size {qc['cut_front_window_size']}",
            f"--cut_front_mean_quality {qc['cut_front_mean_quality']}"
        ])
    if qc["cut_tail"]:
        params.extend([
            f"--cut_tail --cut_tail_window_size {qc['cut_tail_window_size']}",
            f"--cut_tail_mean_quality {qc['cut_tail_mean_quality']}"
        ])
    if qc["cut_right"]:
        params.extend([
            f"--cut_right --cut_right_window_size {qc['cut_right_window_size']}",
            f"--cut_right_mean_quality {qc['cut_right_mean_quality']}"
        ])

    # Length filtering
    lf = config["fastp"]["length_filtering"]
    params.extend([
        f"--length_required {lf['length_required']}",
        f"--length_limit {lf['length_limit']}" if lf['length_limit'] else ""
    ])

    # Overrepresented sequence analysis
    if config["fastp"]["overrepresented_sequence_analysis"]["overrepresentation_analysis"]:
        osa = config["fastp"]["overrepresented_sequence_analysis"]
        params.extend([
            "--overrepresentation_analysis",
            f"--overrepresentation_sampling {osa['overrepresentation_sampling']}"
        ])

    # Concatenate all parameters
    return " ".join(params)

rule trimming:
    input:
        se = 'Preprocessing/{type}.se.fq'
    output:
        se_trimmed = temp('Preprocessing/{type}.se.trimmed.fq'),
        json = 'Stats/{type}/fastp.{type}.json',
        html = 'Stats/{type}/fastp.{type}.html'
    threads:
        getThreads(10)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    params:
        fastp_params = lambda wildcards: get_fastp_params(wildcards)
    log:
        'logs/preprocessing/fastp_trimming.{type}.log'
    message:
        "Fastp trimming: Trimming {wildcards.type} reads."
    shell:
        """
        fastp --in1 {input.se} \
              --out1 {output.se_trimmed} \
              {params.fastp_params} \
              --thread {threads} \
              --json {output.json} \
              --html {output.html} \
              >> {log} 2>&1
        """
