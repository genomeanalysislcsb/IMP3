##################################################
# Mantis
##################################################
# https://github.com/PedroMTQ/mantis
# NOTE: Expecting Mantis to be included as a git submodule
#       Tested with version 14f75ac.

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
    params:
        mod_path=os.path.join(MOD_DIR, "mantis")
    threads: 8
    conda:
        os.path.join(ENV_DIR, "mantis.yaml")
    message:
        "Mantis: Setting up databases and other resources"
    shell:
        "python {params.mod_path}/ setup_databases --mantis_config {input} --cores {threads} &> {log} && "
        "python {params.mod_path}/ check_sql --mantis_config {input} --cores {threads} &>> {log} && "
        "touch {output}"