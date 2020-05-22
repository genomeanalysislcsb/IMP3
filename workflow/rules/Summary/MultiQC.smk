# MultiQC summary report

from os.path import splitext

MQC_MODULES = []

# FastQC
if MG or MT:
    MQC_MODULES.append("fastqc")
    # data type: mg, mt or both
    MQC_FASTQC_T = TYPES
    # read type: PE or SE or both, raw and processed
    MQC_FASTQC_R_RAW  = ["r1", "r2"]
    MQC_FASTQC_R_PROC = ["r1", "r2", "se"]
    if not (MG and MT) and (len(MG) <= 1 or len(MT) <= 1):
        MQC_FASTQC_R_RAW  = ["se"]
        MQC_FASTQC_R_PROC = ["se"]
    else:
        assert len(MG) == 2 and len(MT) == 2
    MQC_FASTQC_R = expand("{rtype}_fastqc", rtype=MQC_FASTQC_R_RAW)
    if 'preprocessing' in IMP_STEPS:
        MQC_FASTQC_R += expand("{rtype}.preprocessed_fastqc", rtype=MQC_FASTQC_R_PROC)
# Prokka
if "analysis" in IMP_STEPS:
    MQC_MODULES.append("prokka")

rule multiqc_report:
    input:
        # FastQC
        fastqc=expand(
            'Stats/{type}/{type}.{rtype}.zip',
            type=MQC_FASTQC_T,
            rtype=MQC_FASTQC_R
        ) if (MG or MT) else [],
        # Prokka
        prokka=["Analysis/annotation/prokka.txt"] if "analysis" in IMP_STEPS else [],
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
        # for file in {input.fastqc} 
        # do
        #     if [ -s ${{file}} ] && [ ! $(unzip -l -q ${{file}} | grep -q \"No errors\") ]
        #     then
        #         echo ${{file}} >> {output.include}
        #     else
        #         echo ${{file}} >> {output.exclude}
        #     fi
        # done
        echo '{input.fastqc} {input.prokka}' | sed 's/\s\+/\\n/g' >> {output.include}
        multiqc --file-list {output.include} {params.modules} --config {params.configfile} --filename {params.filename} &> {log}
        """

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
