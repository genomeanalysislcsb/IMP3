# Utils to be used in Snakemake files

def extend_config(config, config_path):
    """
    Merge given config (obj) with data from another config file (path)
    """
    import os
    if os.path.exists(config_path):
        print("reading configuration from " + config_path)
        with open(config_path, 'r') as rhandle:
            data = yaml.load(rhandle)
            config = dict_merge(config, data)
    return config

def dump_config(config, yaml_path, mode='w+'):
    """
    Dump a YAML object (config) into a file
    """
    import yaml
    from collections import OrderedDict
    yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
    yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))
    with open(yaml_path, mode) as ofile:
        yaml.dump(config, ofile, allow_unicode=True, default_flow_style=False)

def dict_merge(a, b):
    """
    Deep merge 2 dicts together
    """
    from copy import deepcopy
    if not isinstance(b, dict):
        return b
    result = deepcopy(a)
    for k, v in b.items():
        if k in result and isinstance(result[k], dict):
            result[k] = dict_merge(result[k], v)
        else:
            result[k] = deepcopy(v)
    return result

def get_sample(config, outputdir):
    import os
    import re
    sample = os.environ.get("SAMPLE", config['sample'])
    if sample == "":
        sample = "_".join(outputdir.split("/")[-2:])
    sample = re.sub("[;|.-]", "_", sample)
    sample = re.sub("_+", "_", sample)
    return sample

def create_db_dir(dbpath):
    """
    Create DB directory
    """
    import os
    if not os.path.exists(dbpath):
        os.makedirs(dbpath)

def create_tmp_dir(tmpdir, outputdir):
    """
    temporary directory will be stored inside the OUTPUTDIR directory
    unless a absolute path is set
    """
    import os
    if not os.path.isabs(tmpdir):
        tmpdir = os.path.join(outputdir, tmpdir)
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    return tmpdir

def get_hardware(config):
    import os
    bigcores    = os.environ.get("BIGCORENO", config['mem']['big_mem_cores'])
    bigtotal    = os.environ.get("BIGMEMTOTAL", config['mem']['big_mem_total_gb'])
    mem     = str(os.environ.get("MEMCORE", config['mem']['normal_mem_per_core_gb'])) + "G"
    bigmem  = str(os.environ.get("BIGMEMCORE", config['mem']['big_mem_per_core_gb'])) + "G"
    print("BIGMEM: cores: {} w/ {}; total {}\nMEM: {}".format(bigcores, bigmem, bigtotal, mem))
    return bigcores, bigtotal, mem, bigmem

def prepare_input_files(inputs, outputs):
    """
    Prepare file names from input into snakemake pipeline.
    """
    import os
    if len(inputs) != len(outputs):
        raise OSError("//Inputs and outputs are not of the same length: %s <> %s" % (', '.join(inputs), ', '.join(outputs)))
    for infilename, outfilename in zip(inputs, outputs):
        _, fname1 = os.path.split(infilename)
        _process_file(fname1, infilename, outfilename)

def _process_file(fname, inp, outfilename):
    """
    Write the input to the output. Handle raw, zip, or bzip input files.
    """
    print(inp, '=>', outfilename)
    import os
    import bz2
    import shutil
    # ungunzip
    if os.path.splitext(fname)[-1] in ['.gz', '.gzip']:
        with open(outfilename, 'wb') as whandle, gzip.open(inp, 'rb') as rhandle:
            whandle.write(rhandle.read())
    # unbzip2
    elif os.path.splitext(fname)[-1] in ['.bz2', '.bzip2']:
        shell("bzip2 -dc {i} > {o}".format(i=inp, o=outfilename))
    # copy
    else:
        shutil.copy(inp, outfilename)
