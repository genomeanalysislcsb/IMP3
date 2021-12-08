#!/usr/bin/env python3

# NOTE: This script uses some snakemake keywords and expects certain names arguments:
#       snakemake.output.download: directory, required: will be used to save dowloaded files
#       snakemake.output.result: file/directory, optional: will be used copy downloaded and decompressed file(s)
#       snakemake.params.url: URL (str), required: URL to be used to dowload data
#       snakemake.log: file(s), optional; the first specified file will be used for logging
#       snakemake.threads: optional (snakemake's default: 1)

##################################################
# MODULES
##################################################

import os
import logging

##################################################
# FUNCTIONS
##################################################
def setup_logger(logfile=None):
    """Setup a logger"""
    import logging
    # create logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # formatter
    logform = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logform)
    logger.addHandler(ch)
    # file handler
    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setFormatter(logform)
        logger.addHandler(fh)
    return

def create_temp(temp_type, prefix="imp3_install_", dir_path=snakemake.output.download):
    """
    Create a temporary file or directory: "<dir_path>/<prefix>_.*".
    Note: the created file/directory need to be removed by the user.

    Parameters
    ----------
    temp_type : str
        Type: "file" or "directory"
    prefix : str
        Name prefix, set to None if not required
    dir_path : str
        Output path, set to None to use defaul temp path
    
    Returns
    -------
    str or tuple(file handle, str)
        Path for directories, tuple of file handle and path for files
    
    Raises
    ------
    Exception
        If given type is not known, i.e. not "file" or "directory"
    """
    import tempfile
    if temp_type == "file":
        return tempfile.mkstemp(prefix=prefix, dir=dir_path)
    if temp_type == "directory":
        return tempfile.mkdtemp(prefix=prefix, dir=dir_path)
    else:
        raise Exception("Unknown temporary output type: {}".format(temp_type))

def run_cmd(cmd, stdout=None):
    """
    Run given command.

    Parameters
    ----------
    cmd : str
        Command to be executed
    stdout : file handle
        File handle to be used to save stdout; set to None if not required

    Raises
    ------
    Assertion error
        If the return code of the executed command is not 0.
    """
    import shlex
    import subprocess
    logging.info("RUN: {}, stdout={}".format(cmd, stdout))
    cmd_run = subprocess.run(shlex.split(cmd), stdout=stdout, stderr=subprocess.PIPE)
    assert cmd_run.returncode == 0, "Error when running CMD: {}\n{}".format(cmd, cmd_run.stderr)

def get_file_ext(fname):
    """
    Wrapper for filetype.guess()

    Parameters
    ----------
    fname : str
        File path which should be decompressed
    
    Returns
    -------
    obj
        Output of filetype.guess(); None for directories and if no extension could be guessed.
    """
    import filetype
    try:
        fname_ext = filetype.guess(fname)
    except IsADirectoryError: # path is a directory, nothing to decompress
        return None
    return fname_ext

def decompress(fname, threads=snakemake.threads):
    """
    For given path, guess the file type and perform decompression if necessary and implemented.
    If the path is a directory then nothing will be done.
    
    Parameters
    ----------
    fname : str
        File path which should be decompressed
    
    Returns
    -------
    str
        File/directory path after decompression (same as input file path if no decompression was done)
    
    Raises
    ------
    Exception
        If there is no implemented case for the guessed file type extension.
    """
    # file extension
    fname_ext = get_file_ext(fname)
    # no extension, return original name (i.e. no decompression was done)
    if not fname_ext:
        logging.info("No decompression required: {}".format(fname))
        return fname
    else:
        logging.info("Decomporess: {}".format(fname))
    # decompress depending on extension
    if fname_ext.extension == "gz":
        fname_new_handle, fname_new = create_temp(temp_type="file")
        run_cmd(cmd="pigz -p {p} -d -c {i}".format(i=fname, p=threads), stdout=fname_new_handle)
    elif fname_ext.extension == "zip":
        fname_new = create_temp(temp_type="directory")
        run_cmd(cmd="unzip {i} -d {d}".format(i=fname, d=fname_new))
    elif fname_ext.extension == "tar":
        fname_new = create_temp(temp_type="directory")
        run_cmd(cmd="tar -xf {i} -C {d}".format(i=fname, d=fname_new))
    elif fname_ext.extension == "bz2": # NOTE: not tested yet
        fname_new_handle, fname_new = create_temp(temp_type="file")
        run_cmd(cmd="pbzip2 -p {p} -d -c {i}".format(i=fname, p=threads), stdout=fname_new_handle)
    # TODO: other extensions
    else:
        raise Exception("Unknown file extension: {}".format(fname_ext.extension))
    fname_new_ext = get_file_ext(fname_new)
    # check:
    assert fname_new_ext is None or fname_ext.extension != fname_new_ext.extension, "Extension has not changed after decompression: {}".format(fname_new_ext.extension)
    return fname_new

def decompress_recursively(fname, iters_max=4):
    """
    For given path, decompression if necessary as long as decompressin can be applied.
    
    Parameters
    ----------
    fname : str
        File path which should be decompressed
    iters_max : int
        Maximal number of decompression iterations to perform (to avoid an infinte loop if something goes wrong)
    
    Returns
    -------
    str
        File/directory path after complete decompression (same as input file path if no decomporessin was done)
    
    Raises
    ------
    Exception
        If there is no implemented case for the guessed file type extension.
    Assertion error
        If decompression has been performed more than N times
    """
    # to avoid an infinite loop, decompression will be only performed that many times
    # if the number is exceeded, an assertion error will be thrown.
    fname_orig = fname
    iters = 0
    # file name after decompression (same as input path if no decompression is required)
    fname_new = None
    # decompress as long as returned file name is not the same as given file name
    while fname != fname_new:
        if fname_new:
            fname = fname_new
        fname_new = decompress(fname)
        iters += 1
        assert iters < iters_max, "Decompression has been already beformed {} times on the same input file {}. Something is wrong!".format(iters, fname_orig)
    return fname_new

def download(url, target):
    """Download file from given URL into given output file using wget"""
    logging.info("Download: {} -> {}".format(url, target))
    run_cmd(cmd="wget -O {output} {url}".format(url=url, output=target))

def copy_data(source, target):
    """Copy data from source (path) to target (path) using rsync"""
    import os
    # add trailing slash so only the content is copied but the parent folder itself
    if os.path.isdir(source):
        source = os.path.join(source, "")
    logging.info("Copy data: {} -> {}".format(tmp_download, target))
    run_cmd("rsync -rtP {} {}".format(source, target))

##################################################
# MAIN
##################################################
# required snakemake attributes
assert hasattr(snakemake.output, "download"), "There is no \"download\" keyword in output of the snakemake"
assert hasattr(snakemake.params, "url"), "There is no \"url\" keword in params of the snakemake rule"

# logger
if hasattr(snakemake, "log") and len(snakemake.log) > 0:
    setup_logger(logfile=snakemake.log[0])
else:
    setup_logger(logfile=None)

# create folder for downloads
if not os.path.exists(snakemake.output.download):
    os.makedirs(snakemake.output.download)

# temp file for the download
tmp_download = create_temp(temp_type="file")[1]

# download into tmp dir.
download(snakemake.params.url, tmp_download)

# decompress in tmp dir.
tmp_decompressed = decompress_recursively(tmp_download)

# create output
if hasattr(snakemake.output, "result"):
    copy_data(tmp_decompressed, snakemake.output.result)