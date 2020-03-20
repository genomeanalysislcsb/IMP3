# include rules for the Report part of the workflow

include:
    "../../../Report/diagramm.smk"

include:
    "../../../Report/libraries.smk"

include:
    "../../../Report/report.smk"

# master command
rule REPORT:
    input:
        "REPORT.html"
    output:
        touch('report.done')
