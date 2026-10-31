
arch_inputs = "status/workflow.done"
if MG or MT:
    arch_input.append("status/reads.compressed")

localrules: ARCHIVING

rule ARCHIVING:
    input:
        arch_input
    output:
        "status/archiving.done"

rule sync_final:
    input:
        OUTPUTDIR,
        "status/workflow.done"
    output:
        SUMDIR
    threads:
        1
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    message:
        "Synching {OUTPUTDIR} to {SUMDIR}"
    shadow: "copy-minimal"
    shell:
        """
        rsync -a \
              --remove-source-files \
              {input[0]} \
              {output}
        """

def entry_gz:
    if "preprocessing" in IMP_STEPS:
        ctrl.append("status/preprocessing.done")
    elif "summary" in SUM_STEPS:
        ctrl.append("status/stats.done")
    if "assembly" in IMP_STEPS and (not MGaln and not MTaln):
        ctrl.append("status/assembly.done")
    if "taxonomy" in IMP_STEPS:
        ctrl.append("status/taxonomy.done")
    return ctrl
    
def ctrl_readgz(wildcards):
    checkpoint_output = checkpoints.workflow.get(**wildcards).output[0]
    pot_reads = expand("Preprocessing/{i}.fq.gz", i=glob_wildcards("Preprocessing/{i}.fq").i)
    real_reads = list(set(pot_reads).difference([expand("Preprocessing/{type}.{part}.preprocessed.fq.gz",
                  type=TYPES,part=["r1","r2","se"]),
                                                 expand("Preprocessing/{type}.{part}.fq.gz",
                  type=TYPES,part=["r1","r2"])]))
    return real_reads

localrules: compress_reads

rules compress_reads:
    input: 
        ctrl_readgz
    output:
        touch("status/reads.compressed")

if MG or MT:
    rule gzip_reads:
        input:
            # Make sure its done at the very end
            check = "status/workflow.done",
            file = "{file}",
            checks = entry_gz  # TODO: ???
        output:
            "{file}.gz"
        threads:
            1
        resources:
            runtime = SHORTISH_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        message:
            "gzip: gzipping {wildcards.file}."
        shadow: "copy-minimal"
        shell:
            "gzip {input.file}"

rule tar_gzip_directory:
    input:
        check="status/workflow.done",
        directory="{directory}"
    output:
        "{directory}.tar.gz"
    threads:
        1
    resources:
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    message:
        "tar_gzip: tarring and gzipping {wildcards.directory}."
    shadow: "copy-minimal"
    shell:
        "tar -cvzf {output} {input.directory} && rm -R {input.directory}"
