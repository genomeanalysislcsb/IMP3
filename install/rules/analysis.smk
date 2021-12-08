##################################################
# Mantis
##################################################
# https://github.com/PedroMTQ/mantis
# NOTE: Expecting Mantis to be included as a git submodule
#       Tested with version 14f75ac.

localrules: mantis_config
rule mantis_config:
    output:
        "mantis/mantis.config"
    params:
        config=config["mantis"]
    message:
        "Mantis: create config file"
    run:
        with open(output[0], "w+") as ofile:
            # default HMMs
            for ref_name, ref_path in params.config["default"].items():
                ofile.write(f'{ref_name}={ref_path}\n')
            # custom HMMs
            if params.config["custom"]:
                for ref_path in params.config["custom"]:
                    ofile.write(f'custom_ref={ref_path}\n')
            # weights
            for weights_name, weights_value in params.config["weights"].items():
                ofile.write(f'{weights_name}={weights_value}\n')

rule mantis_setup:
    input:
        "mantis/mantis.config"
    output:
        "mantis/setup.done"
    log:
        "mantis/setup.log"
    params:
        mod_path=os.path.join(MOD_DIR, "mantis")
    threads: 12
    conda:
        os.path.join(ENV_DIR, "mantis.yaml")
    message:
        "Mantis: Setting up databases and other resources"
    shell:
        "python {params.mod_path}/ setup_databases --mantis_config {input} --cores {threads} &> {log} && touch {output}"