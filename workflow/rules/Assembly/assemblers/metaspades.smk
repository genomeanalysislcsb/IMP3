# metaSPAdes
#   manual: https://cab.spbu.ru/files/release3.13.0/manual.html

# TODO: LR indep. of the omic?

rule metaspades:
    input:
        r1="Preprocessing/{omic}.r1.preprocessed.fq",
        r2="Preprocessing/{omic}.r2.preprocessed.fq",
        se="Preprocessing/{omic}.se.preprocessed.fq",
        lr=metaspades_lr_input # TODO: indep. of omic ???
    output:
        fa="Assembly/intermediary/{omic}.metaspades_preprocessed.1/contigs.fasta",
        fa_link="Assembly/intermediary/{omic}.metaspades_preprocessed.1.fa"
    log:
        "logs/assembly_metaspades_assembly_from_preprocessing.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics()),
    params:
        odir=lambda wildcards, output: os.path.dirname(output.fa),
        mem=get_mem_big(get_cores_big(), "Gb"), # mem in Gb
        kmers=",".join(config["assembly"]["spades"]["metaspades"]["kmers"]),
        lr=metaspades_lr_input_params
    threads: get_cores_big()
    resources:
        runtime="120:00:00",
        mem=get_mem_big(get_cores_big()),
    conda:
        os.path.join(ENVDIR, "IMP_assembly.yaml")
    message:
        "MetaSPAdes assembly: {input}"
    shell:
        # clean up if needed
        "if [ -d {params.odir} ]; then rm -rf {params.odir}; fi && "
        # assembly
        "spades.py --meta --pe1-1 {input.r1} --pe1-2 {input.r2} --pe1-s {input.se} {params.lr} {input.lr} -o {params.odir} "
        "-t {threads} -m {params.mem} -k {params.kmers} &> {log} && "
        # link to output
        "ln -fsv $(realpath --relative-to=\"$(dirname {output.fa_link})\" {output.fa}) {output.fa_link} && "
        "touch -h {output.fa_link}" # TODO: why ???

use rule metaspades as metaspades_unmapped with:
    input:
        r1="Assembly/intermediary/{omic}.r1.unmapped.fq",
        r2="Assembly/intermediary/{omic}.r2.unmapped.fq",
        se="Assembly/intermediary/{omic}.se.unmapped.fq",
        lr=[],
    output:
        fa="Assembly/intermediary/{omic}.metaspades_unmapped.2/contigs.fasta",
        fa_link="Assembly/intermediary/{omic}.metaspades_unmapped.2.fa",
    log:
        "logs/assembly_metaspades_assembly_from_unmapped.{omic}.log"

# TODO: https://git-r3lab.uni.lu/IMP/imp3/-/issues/36