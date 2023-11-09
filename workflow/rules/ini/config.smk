# settings
# Importing steps from the configuration file and sorting them
IMP_STEPS = sorted(config['steps'].split())
# Importing preprocessing_filtering from the configuration file
PREPROCESSING_FILTERING = config['preprocessing_filtering']
# Importing filter from the configuration file and splitting it
FILTER = config['filtering']['filter'].split()
# Importing assembler from the configuration file
IMP_ASSEMBLER = config['assembly']['assembler']
# Importing binners from the configuration file and splitting them
BIN_STEPS = config['binning']['binners'].split()
# Checking if the number of binning methods exceeds the maximum limit of 6
if len(BIN_STEPS) > 6:
    raise Exception('A maximum of six binning methods per run are supported at the moment.')
# If neither MG nor MGaln is set, set BIN_STEPS to ['MetaBAT']
if not MG and not MGaln:
    BIN_STEPS = ['MetaBAT']
# Importing summary_steps from the configuration file and splitting them
SUM_STEPS = config['summary_steps'].split()

# hardware parameters
# Defining dictionaries for byte conversion
byte_dict = {'K': 1e3, 'M': 1e6, 'G': 1e9, 'T': 1e12}
to_gb_dict = {'K': 1e-6, 'M': 1e-3, 'G': 1e0, 'T': 1e3}
to_mb_dict = {'K': 1e-3, 'M': 1e0, 'G': 1e3, 'T': 1e6}

# Importing cluster specifications from the configuration file
MAIN_BIGCORENO = config['cluster_specs']['main_cluster']['big_mem_cores']
MAIN_BIGMEMTOTAL = config['cluster_specs']['main_cluster']['big_mem_total']
MAIN_MEMCORE = str(config['cluster_specs']['main_cluster']['mem_per_core']) + config['cluster_specs']['mem_unit']
MAIN_BIGMEMCORE = str(config['cluster_specs']['main_cluster']['big_mem_per_core']) + config['cluster_specs']['mem_unit']
MAIN_CORESMAX = config['cluster_specs']['main_cluster']['cores_max']
MAIN_BIGCORESMAX = config['cluster_specs']['main_cluster']['big_mem_cores_max']

# Converting memory units to MB
MAIN_MEMCORE_MB = int(float(MAIN_MEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])
MAIN_BIGMEMCORE_MB = int(float(MAIN_BIGMEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])

# Checking if high_mem_cluster is set in the configuration file
if config['cluster_specs'].get('high_mem_cluster'):
    print('#' * 10, 'HIGH MEM ON', '#' * 10)
    USE_HIGH_MEM_CLUSTER = True
    # Importing high_mem_cluster specifications from the configuration file
    HIGH_MEM_BIGCORENO = config['cluster_specs']['high_mem_cluster']['big_mem_cores']
    HIGH_MEM_BIGMEMTOTAL = config['cluster_specs']['high_mem_cluster']['big_mem_total']
    HIGH_MEM_MEMCORE = str(config['cluster_specs']['high_mem_cluster']['mem_per_core']) + config['cluster_specs']['mem_unit']
    HIGH_MEM_BIGMEMCORE = str(config['cluster_specs']['high_mem_cluster']['big_mem_per_core']) + config['cluster_specs']['mem_unit']
    HIGH_MEM_CORESMAX = config['cluster_specs']['high_mem_cluster']['cores_max']
    HIGH_MEM_BIGCORESMAX = config['cluster_specs']['high_mem_cluster']['big_mem_cores_max']

    # Converting high_mem_cluster memory units to MB
    HIGH_MEM_MEMCORE_MB = int(float(HIGH_MEM_MEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])
    HIGH_MEM_BIGMEMCORE_MB = int(float(HIGH_MEM_BIGMEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])
else:
    USE_HIGH_MEM_CLUSTER = False

# Checking if the temporary directory is an absolute path
TMPDIR = config['tmp_dir']
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.join(OUTPUTDIR, TMPDIR)
# Creating the temporary directory if it does not exist
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

# annotation
# If annotation is set to 'mantis' in the configuration file, set hmm_DBs to an empty string
if config['annotation'] == 'mantis':
    config['hmm_DBs'] = ''