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
        def write_ref_path(file_handle, ref_path, ref_name):
            """Write ref path to file"""
            if ref_path:
                assert os.path.exists(ref_path) and os.path.isdir(ref_path), f"Mantis ref path does not exist or is not a path: {ref_path}"
                file_handle.write(f"{ref_name}={ref_path}\n")
            elif ref_name != "custom_ref":
                file_handle.write(f"{ref_name}=NA\n")
        
        with open(output[0], "w+") as ofile:
            # default HMMs
            for ref_name, ref_path in params.config["default"].items():
                write_ref_path(ofile, ref_path, ref_name)
            # custom HMMs
            if params.config["custom"]:
                for ref_path in params.config["custom"]:
                    write_ref_path(ofile, ref_path, "custom_ref")
            # weights
            for weights_name, weights_value in params.config["weights"].items():
                ofile.write(f"{weights_name}={weights_value}\n")

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