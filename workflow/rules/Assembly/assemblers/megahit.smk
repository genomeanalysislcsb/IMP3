rule megahit:
    input:
        r1="Preprocessing/{omic}.r1.preprocessed.fq",
        r2="Preprocessing/{omic}.r2.preprocessed.fq",
        se="Preprocessing/{omic}.se.preprocessed.fq",
    output:
        fa="Assembly/intermediary/{omic}.megahit_preprocessed.1/final.contigs.fa",
        fa_link="Assembly/intermediary/{omic}.megahit_preprocessed.1.fa",
    log:
        "logs/assembly_megahit_assembly_from_preprocessing.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics()),
    params:
        odir=lambda wildcards, output: os.path.dirname(output.fa),
        mem=lambda wildcards, threads: get_mem_big(wildcards, None, threads, None, "B"), # mem in bytes
        kmin=config["assembly"]["megahit"]["kmin"],
        kmax=config["assembly"]["megahit"]["kmax"],
        kstep=config["assembly"]["megahit"]["kstep"],
    threads: get_cores_big()
    resources:
        runtime="120:00:00",
        mem=get_mem_big,
    conda:
        os.path.join(ENVDIR, "IMP_assembly.yaml")
    message:
        "MEGAHIT assembly: {input}"
    shell:
        # clean up if needed
        "if [ -d {params.odir} ]; then rm -rf {params.odir}; fi && "
        # assembly
        "megahit -1 {input.r1} -2 {input.r2} -r {input.se} -o {params.odir} "
        "--k-min {params.kmin} --k-max {params.kmax} --k-step {params.kstep} --bubble-level 0 " # TODO: why "--bubble-level 0"?
        "-t {threads} --cpu-only -m {params.mem} --mem-flag 1 &> {log} && "
        # link to output
        "ln -fsv $(realpath --relative-to=\"$(dirname {output.fa_link})\" {output.fa}) {output.fa_link} && "
        "touch -h {output.fa_link}" # TODO: why ???

use rule megahit as megahit_unmapped with:
    input:
        r1="Assembly/intermediary/{omic}.r1.unmapped.fq",
        r2="Assembly/intermediary/{omic}.r2.unmapped.fq",
        se="Assembly/intermediary/{omic}.se.unmapped.fq",
    output:
        fa="Assembly/intermediary/{omic}.megahit_unmapped.2/final.contigs.fa",
        fa_link="Assembly/intermediary/{omic}.megahit_unmapped.2.fa"
    log:
        "logs/assembly_megahit_assembly_from_unmapped.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics())

rule megahit_mgmt:
    input:
        r1=expand("Preprocessing/{omic}.r1.preprocessed.fq", omic=["mg", "mt"]),
        r2=expand("Preprocessing/{omic}.r2.preprocessed.fq", omic=["mg", "mt"]),
        se=expand("Preprocessing/{omic}.se.preprocessed.fq", omic=["mg", "mt"]) + [
            # MT assembly (primary + from unmapped reads)
            "Assembly/intermediary/mt.megahit_preprocessed.1/final.contigs.fa", # TODO: why not the link?
            "Assembly/intermediary/mt.megahit_unmapped.2/final.contigs.fa", # TODO: why not the link?
        ],
    output:
        fa="Assembly/intermediary/mgmt.megahit_hybrid.1/final.contigs.fa",
        fa_link="Assembly/intermediary/mgmt.megahit_hybrid.1.fa",
    log:
        "logs/assembly_megahit_hybrid_assembly_1.log"
    params:
        odir=lambda wildcards, output: os.path.dirname(output.fa),
        mem=lambda wildcards, threads: get_mem_big(wildcards, None, threads, None, "B"), # mem in bytes
        kmin=config["assembly"]["megahit"]["kmin"],
        kmax=config["assembly"]["megahit"]["kmax"],
        kstep=config["assembly"]["megahit"]["kstep"],
    threads: get_cores_big()
    resources:
        runtime="120:00:00",
        mem=get_mem_big,
    conda:
        os.path.join(ENVDIR, "IMP_assembly.yaml")
    message:
        "MEGAHIT: {input}"
    message:
        "MEGAHIT MG/MT assembly: {input}"
    shell:
        # clean up if needed
        "if [ -d {params.odir} ]; then rm -rf {params.odir}; fi && "
        # assembly
        "megahit -1 $(echo {input.r1} | tr ' ' ',') -2 $(echo {input.r2} | tr ' ' ',') -r $(echo {input.se} | tr ' ' ',') -o {params.odir} "
        "--k-min {params.kmin} --k-max {params.kmax} --k-step {params.kstep} "
        "-t {threads} --cpu-only -m {params.mem} --mem-flag 1 &> {log} && "
        # link to output
        "ln -fsv $(realpath --relative-to=\"$(dirname {output.fa_link})\" {output.fa}) {output.fa_link} && "
        "touch -h {output.fa_link}" # TODO: why ???

use rule megahit_mgmt as megahit_mgmt_unmapped with:
    input:
        r1=expand("Assembly/intermediary/mgmt.r1.megahit_hybrid.{omic}.unmapped.fq", omic=["mg", "mt"]),
        r2=expand("Assembly/intermediary/mgmt.r2.megahit_hybrid.{omic}.unmapped.fq", omic=["mg", "mt"]),
        se=expand("Assembly/intermediary/mgmt.se.megahit_hybrid.{omic}.unmapped.fq", omic=["mg", "mt"]),
    output:
        fa="Assembly/intermediary/mgmt.megahit_hybrid.2/final.contigs.fa",
        fa_link="Assembly/intermediary/mgmt.megahit_hybrid.2.fa",
    log:
        "logs/assembly_megahit_hybrid_assembly_2.log"
