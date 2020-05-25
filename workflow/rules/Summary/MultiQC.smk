# MultiQC summary report

from os.path import splitext, isfile

##############################
# Help variables
MQC_MODULES = []
# FastQC
if MG or MT:
    MQC_MODULES.append("fastqc")
    # read type: PE or SE or both, raw and processed
    MQC_FASTQC_RTYPES_RAW  = ["r1", "r2"]
    MQC_FASTQC_RTYPES_PROC = ["r1", "r2", "se"]
    if not (MG and MT) and (len(MG) <= 1 or len(MT) <= 1):
        MQC_FASTQC_RTYPES_RAW  = ["se"]
        MQC_FASTQC_RTYPES_PROC = ["se"]
    else:
        assert len(MG) == 2 and len(MT) == 2
    MQC_FASTQC_RTYPES = expand("{rtype}_fastqc", rtype=MQC_FASTQC_RTYPES_RAW)
    if "preprocessing" in IMP_STEPS:
        MQC_FASTQC_RTYPES += expand("{rtype}.preprocessed_fastqc", rtype=MQC_FASTQC_RTYPES_PROC)
# Trimmomatic/CutAdapt
if 'preprocessing' in IMP_STEPS and (not config['nextseq'] or not isfile(DBPATH + "/adapters/adapters.done")):
    # CutAdapt: TODO
    if config['nextseq'] and not isfile(DBPATH + "/adapters/adapters.done"):
        MQC_MODULES.append("cutadapt")
        # MQC_CUTADAPT_RTYPES = [] # TODO
    # Trimmomatic
    else:
        MQC_MODULES.append("trimmomatic")
# SortMeRNA
if PREPROCESSING_FILTERING and MT and len(MT) > 1:
    MQC_MODULES.append("sortmerna")
# Prokka
if "analysis" in IMP_STEPS:
    MQC_MODULES.append("prokka")
# featureCounts
if "analysis" in IMP_STEPS:
    MQC_MODULES.append("featureCounts")
    MQC_FC_CTYPES = ["%s_counts" % s for s in ["rRNA","CDS"]] + ["%s.counts" % s for s in config["hmm_DBs"].split()]

##############################
# Input processing
localrules: multiqc_freaturecounts

rule multiqc_freaturecounts:
    input:
        "Analysis/annotation/{type}.{ctype}.tsv.summary"
    output:
        temp("Analysis/annotation/{type}.{ctype}.multiqc.summary")
    wildcard_constraints:
        type="|".join(TYPES),
        ctype="|".join(MQC_FC_CTYPES)
    shell:
        """
        ctype=$(echo \"{wildcards.ctype}\" | sed \"s/\(\.\|_\)counts//\");
        sed \"s/^Status\\t.*/Status\\t{wildcards.type}.${{ctype}}/\" {input} > {output}
        """

##############################
# MultiQC report
rule multiqc_report:
    input:
        # FastQC
        fastqc=expand(
            'Stats/{type}/{type}.{rtype}.zip',
            type=TYPES,
            rtype=MQC_FASTQC_RTYPES
        ) if "fastqc" in MQC_MODULES else [],
        # CutAdapt: TODO
        # cutadapt=expand(
        #     "logs/preprocessing_nextseq_trimming.{type}.{rtype}.log",
        #     type=TYPES,
        #     rtype=MQC_CUTADAPT_RTYPES
        # ) if "cutadapt" in MQC_MODULES else [],
        # Trimmomatic
        trimmomatic=expand(
            "logs/preprocessing_trimming.{type}.log",
            type=TYPES
        ) if "trimmomatic" in MQC_MODULES else [],
        # SortMeRNA
        sortmerna=[
            "Preprocessing/mt.pe.trimmed.rna.log",
            "Preprocessing/mt.se.trimmed.rna.log"
        ] if "sortmerna" in MQC_MODULES else [],
        # Prokka
        prokka=[
            "Analysis/annotation/prokka.txt"
        ] if "prokka" in MQC_MODULES else [],
        # featureCounts
        featurecounts=expand(
            "Analysis/annotation/{type}.{ctype}.multiqc.summary",
            type=TYPES,
            ctype=MQC_FC_CTYPES
        ) if "featureCounts" in MQC_MODULES else [],
    output:
        include='Stats/multiqc_files.txt',
        exclude='Stats/multiqc_files_excluded.txt',
        html='Stats/multiqc_report.html',
        data='Stats/multiqc_report_data.zip'
    log:
        "logs/multiqc.log"
    threads: 1
    params:
        runtime = "00:30:00",
        mem = MEMCORE,
        configfile = CONFIGDIR + "/config.multiqc.yaml",
        modules = " ".join(["--module %s" % m for m in MQC_MODULES]),
        filename = lambda wildcards, output: splitext(output.html)[0]
    conda:
        ENVDIR + "/MultiQC.yaml"
    message:
        "Create MultiQC report from {input}"
    shell:
        """
        > {output.include}
        > {output.exclude}
        echo '{input}' | sed 's/\s\+/\\n/g' >> {output.include}
        multiqc --file-list {output.include} {params.modules} --config {params.configfile} --filename {params.filename} &> {log}
        """

# MultiQC report status
localrules: multiqc_done
rule multiqc_done:
    input:
        html='Stats/multiqc_report.html',
        data='Stats/multiqc_report_data.zip'
    output:
        "status/multiqc.done"
    threads: 1
    message:
        "MultiQC report has been created: {output}"
    shell:
        "touch {output}"
