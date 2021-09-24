##################################################
# CONFIG VALIDATION

# default config file: any other config file or config parameter from the CMD will overwride its parameters
configfile:
    srcdir("../../../config/config.yaml")

# validate the config file
validate(config, schema="../../../schemas/config.schema.yaml")

if config["input"]["lr"]["mg"]:
    assert lr_tech, "Long-read data given bot long-read seq. technology is not set."

if not config["sample"]:
    config["sample"] = get_sample_name(OUTPUTDIR)

if not has_sr("mg") and not has_asm_map("mg") and "binning" in config["steps"]:
    assert "MaxBin" not in config["binning"]["binners"], "MaxBin should be used but there is no MG reads or MG BAM files"
    assert "binny" not in config["binning"]["binners"], "binniny should be used but there is no MG reads or MG BAM files"

if has_mgmt():
    assert has_sr("mg", "pe") == has_sr("mt", "pe"), "MG and MT reads must have the same read types: PE and/or SE"
    assert has_sr("mg", "se") == has_sr("mt", "se"), "MG and MT reads must have the same read types: PE and/or SE"

##################################################
# PATHS

# paths
# include global functions
ROOTDIR = srcdir(".")
SRCDIR = srcdir("../../scripts")
BINDIR = srcdir("../../bin")
ENVDIR = srcdir("../../envs")
MODDIR = srcdir("../../../submodules")
# output
OUTPUTDIR = get_abs_path(config['out_path'])
# summary
SUMDIR = os.path.join(OUTPUTDIR, "summary")
# databases
DBPATH = get_abs_path(config['db_path'])
if not os.path.exists(DBPATH):
    os.makedirs(DBPATH)
# temporary directory will be stored inside the OUTPUTDIR directory
# unless a absolute path is set
# TMPDIR = get_abs_path(config['tmp_dir'])
# if not os.path.exists(TMPDIR):
#     os.makedirs(TMPDIR)

# EMAIL = config['email']
# if EMAIL != "":
#     if not re.fullmatch(r"^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$", EMAIL):
#         EMAIL = ""
#         print("Your email address is not valid, you will not receive notifications.")


config["input"] = proc_input_files()
check_input_files()


# TEST
print("Input files: {}".format(get_ori_files()))
print("Omics: {}".format(get_omics()))
print("Reads: {}".format(get_read_types_raw("mg")))
print("Assembly omics: {}".format(get_omics_assembly()))
# print_config()
