##################################################
# XXX

# default executable for snakemake
shell.executable("bash")

##################################################
# CONFIG VALIDATION

# YAML representers
add_config_representer()

# default config file
# NOTE: any other config file or config parameter from the CMD will overwride its parameters
configfile:
    srcdir("../../../config/config.yaml")

# validate the config file
validate(config, schema="../../../schemas/config.schema.yaml")

if config["input"]["lr"]["mg"]:
    assert lr_tech, "Long-read data given but sequencing technology (\"lr_tech\") is not set."

if not config["sample"]:
    config["sample"] = get_sample_name(OUTPUTDIR)

if not has_sr("mg") and not has_asm_map("mg") and "binning" in config["steps"]:
    assert "MaxBin" not in config["binning"]["binners"], "MaxBin should be used but no MG reads or BAM was given"
    assert "binny" not in config["binning"]["binners"], "binniny should be used but no MG reads or BAM was given"

if has_mgmt():
    assert has_sr("mg", "pe") == has_sr("mt", "pe"), "MG and MT reads must have the same read types: PE and/or SE"
    assert has_sr("mg", "se") == has_sr("mt", "se"), "MG and MT reads must have the same read types: PE and/or SE"

# TODO: expand

##################################################
# PATHS

# code
ROOTDIR = srcdir(".")
SRCDIR  = srcdir("../../scripts")
BINDIR  = srcdir("../../bin")
ENVDIR  = srcdir("../../envs")
MODDIR  = srcdir("../../../submodules")

# output
OUTPUTDIR = get_abs_path(config['out_path'])

# summary
SUMDIR = os.path.join(OUTPUTDIR, "summary")

# databases
DBPATH = get_abs_path(config['db_path'])

# tmp
# NOTE: tmpdir is set via snakemake's variable "tmpdir" in "default-resources"

# EMAIL = config['email']
# if EMAIL != "":
#     if not re.fullmatch(r"^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$", EMAIL):
#         EMAIL = ""
#         print("Your email address is not valid, you will not receive notifications.")


# input
config["input"] = proc_input_files()
check_input_files()

##################################################
# NEW WORKING DIRECTORY

# working/output directory
workdir:
    OUTPUTDIR

# save config
save_config(os.path.join(OUTPUTDIR, "sample.config.yaml"))
