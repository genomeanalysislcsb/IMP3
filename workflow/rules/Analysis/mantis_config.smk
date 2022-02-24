# NOTE: This file is also included in the installation workflow.
#       Do not add any other code here if this code is not supposed to be used during the installation as well.

# Mantis: create config file from snakemake config
localrules: mantis_config
rule mantis_config:
    output:
        config_file="mantis.config"
    params:
        config=config["mantis"]
    message:
        "Mantis: config file"
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
