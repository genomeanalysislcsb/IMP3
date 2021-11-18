##################################################
# UTILS: functions
##################################################

def flatten_list(l):
    """Flatten list (or another iterable) of iterables and return new list with sorted and unique elements"""
    return sorted(set(list(chain.from_iterable(l))))

def smk_get_genome_url(wildcards):
    """For snakemake: get URL for given genome name"""
    return config["genomes"][wildcards.genome]["url"]

def get_genome_download_target(genome, preped=False):
    """
    Get snakemake output file for given genome

    Parameters
    ----------
    genome : str
        Genome keyword/name from the config attribute "genomes"
    preped: bool
        Whether the output file is expected to be processed, e.g. decompressed
        (then compression extension will not be added)

    Returns
    -------
    str
        File path for the output file based on given name and expected file extension(s)
    """
    ext_c = config["genomes"][genome]["compression"]
    # original file is compressed and the output is not expected to be decompressed
    if ext_c and not preped:
        return "downloads/files/{{genome}}.{{ext}}.{ext_c}".format(ext_c=ext_c)
    # otherwise
    else:
        return "downloads/files/{genome}.{ext}"

def smk_get_genome_download_target(wildcards):
    """For snakemake: get snakemake output file for given genome"""
    return get_genome_download_target(wildcards.genome, False)

def get_kraken2db_file_names():
    """
    Basenames of files required for a kraken2 DB
    Reference: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases
    """
    return ["taxo.k2d", "opts.k2d", "hash.k2d"]

def smk_get_kraken2db_url(wildcards):
    """For snakemake: get URL for given Kraken2 databdase"""
    return config["kraken2"][wildcards.kdb]["url"]

def get_kraken2db_download_target(kdb):
    """
    Get snakemake output file for a Kraken2 database

    Parameters
    ----------
    kdb : str
        Kraken2 database keyword/name from the config attribute "kraken2"
    
    Returns
    -------
    str
        File path for the output file based on given name
    """
    return "downloads/files/{{kdb}}.{ext}".format(ext=config["kraken2"][kdb]["compression"])

def smk_get_kraken2db_download_target(wildcards):
    """For snakemake: get snakemake output file for a Kraken2 database"""
    return get_kraken2db_download_target(wildcards.kdb)