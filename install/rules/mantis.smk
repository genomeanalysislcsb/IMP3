##################################################
# Mantis
##################################################
# https://github.com/PedroMTQ/mantis
# now uses mantis conda install
# tested with mantis version 1.4.5 
#

include:
    "../../workflow/rules/Analysis/mantis_config.smk"

localrules: mantis_config_install
use rule mantis_config as mantis_config_install with:
    output:
        "mantis/mantis.config"

rule mantis_setup:
    input:
        "mantis/mantis.config"
    output:
        "mantis/setup.done"
    log:
        "mantis/setup.log"
    threads: 8
    conda:
        os.path.join(ENV_DIR, "IMP_mantis.yaml")
    message:
        "Mantis: Setting up databases and other resources"
    shell:
        "mantis setup --mantis_config {input} --cores {threads} &> {log} && "
        "mantis check_sql --mantis_config {input} --cores {threads} &>> {log} && "
        "touch {output}"
