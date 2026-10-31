##################################################
# FUNCTIONS
##################################################

def flatten_list(l):
    """Flatten an iterable of iterables and return a sorted list of unique elements"""
    return sorted(set(list(chain.from_iterable(l))))

def smk_get_genome_url(wildcards):
    """For snakemake: get URL for given genome name"""
    return config["genomes"][wildcards.genome]["url"]

def get_kraken2db_file_names():
    """
    Lits of files (without their path) required for a Kraken2 database
    Reference: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases
    """
    return ["taxo.k2d", "opts.k2d", "hash.k2d"]

def smk_get_kraken2db_url(wildcards):
    """For snakemake: get URL for given Kraken2 databdase"""
    return config["kraken2"][wildcards.kdb]["url"]

def smk_get_git_url(wildcards):
    """For snakemake: get URL for given git repository"""
    return config[wildcards.repo]["url"]

def smk_get_git_commit(wildcards):
    """For snakemake: get commit ID for given git repository"""
    return config[wildcards.repo]["commit"]

##################################################
# RULES
##################################################
# Adapters from Trimmomatic repository
rule git_download:
    output:
        repo=temp(directory("downloads/git/{repo}")),
    params:
        url=smk_get_git_url,
        commit=smk_get_git_commit,
    message:
        "Git: download {params.url}, commit {params.commit}"
    shadow: "copy-minimal"
    shell:
        "git clone {params.url} {output} && "
        "cd {output} && "
        "git checkout {params.commit}"
