#!/usr/bin/env python3

##############################
# IMPORTS
##############################
import pandas
import logging
from timeit import default_timer as timer

##############################
# FUNCTIONS
##############################

def setup_logger(logfile=None, level=logging.DEBUG):
    """
    Setup a logger
    """

    # Formatting
    formatter = logging.Formatter('%(asctime)s - %(message)s', '%d/%m/%Y %I:%M:%S %p')

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(level)

    # Console handler
    ch = logging.StreamHandler()
    logger.addHandler(ch)

    # File handler
    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return


def read_in_gff(fname):
    """Read in a GFF"""
    # Col names (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#description-of-the-format)
    gff_cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr"]
    # read in
    gff = pandas.read_csv(fname, sep="\t", header=None, names=gff_cols, comment="#")
    # NOTE: Any GFF feature with an entry in the FAA file (used for annotation w/ Mantis) is expected to have an ID.
    #       Any entry with an ID is removed as the IDs are used as index/row names.
    #       Cases w/o an ID are reported as this should ideally not happen.
    gff = gff.assign(fID=gff.attr.str.extract(r"(^|;)ID=(?P<id>\w+)(;|$)")["id"])
    if any(pandas.isnull(gff.fID)):
        logging.warning(f"Some features in {fname} don't have an ID:\n{gff.loc[pandas.isnull(gff.fID)]}")
        gff = gff.loc[pandas.notnull(gff.fID), :]
    # set row names: feature ID extracted from the attributes column
    gff.set_index(keys="fID", drop=True, inplace=True, verify_integrity=True)
    return gff


def read_mantis_consensus(fname):
    """Read Mantis consensus table"""
    # read in
    annot = pandas.read_csv(fname, sep="\t", header=0, index_col=0)
    return annot


def mantis_links2df(links_str, sep="|"):
    """
    Split Mantis' protein annotations and return a table with keys and values

    Parameters
    ----------
    links_str : str
        Links string provided by Mantis in the consensus output table
    sep : str
        Character/string used to separate anotations

    Returns
    -------
    pandas.DataFrame
        Table with all annotations (except description) stored in columns 'Key' and 'Value'
    """
    return pandas.DataFrame([tuple(hit.split(":")) for hit in links_str.split(sep)], columns=["Key", "Value"])


def mantis_links_generator(links_str, sep="|"):
    """
    Generator yielding summarized Mantis' annotations in GFF attribute format

    Parameters
    ----------
    links_str : str
        Links string provided by Mantis in the consensus output table
    sep : str
        Character/string used to separate annotations

    Yields
    ------
    str
        Summarized annotation per type (key value) in GFF attribute format (<key>=<value1,value2,...>)
    """
    # table of keys and values w/o description
    df_links = mantis_links2df(links_str, sep)
    df_links = df_links.loc[df_links.Key != "description",:]
    # generator
    for gr_name, gr_df in df_links.groupby("Key"):
        yield "{}={}".format(gr_name, ",".join(gr_df.Value.astype(str)))


def mantis_links_pairs(links_str, sep="|"):
    """
    Mantis' annotations for a protein as a table with key and key/value string in GFF format

    Parameters
    ----------
    links_str : str
        Links string provided by Mantis in the consensus output table
    sep : str
        Character/string used to separate anotations

    Yields
    ------
    pandas
        Annotations type (key value) and key/value pairs in GFF attribute format (<key>=<value>)
    """
    # table of keys and values w/o description
    df_links = mantis_links2df(links_str, sep)
    df_links = df_links.loc[df_links.Key != "description", :]
    # create pairs and drop values
    df_links["Pair"] = df_links.apply("=".join, axis=1)
    df_links.drop("Value", inplace=True, axis=1)
    return df_links

##############################
# MAIN
##############################

if __name__ == '__main__':
    log = snakemake.log[0]
    setup_logger(log, level=logging.INFO)

    # read in
    logging.info('Reading Prokka gff.')
    gff = read_in_gff(snakemake.input.gff)
    logging.info('Reading Mantis output.')
    mantis = read_mantis_consensus(snakemake.input.mantis)

    # check attribute keys appearing in given input files: there should be no overlap
    logging.info('Assure no overlap in input file keys.')
    gff_keys = gff.attr.apply(lambda x: set([y.split("=")[0] for y in x.split(";")])).values
    gff_keys = set.union(*gff_keys)
    mantis_keys = mantis.Links.apply(lambda x: set([y.split(":")[0] for y in x.split(snakemake.params.sep)])).values
    mantis_keys = set.union(*mantis_keys)
    assert len(gff_keys.intersection(mantis_keys)) == 0,\
           f'There are common attribute keys in {input}: {gff_keys.intersection(mantis_keys)}'

    # add GFF attributes based on Mantis' annotations
    logging.info('Add GFF attributes based on Mantis annotations and GFF attributes per db for featureCounts.')

    pids_processed = 0
    times_per_gff_pid = []
    times_per_mantis_pid_total = []
    times_per_mantis_pid_links = []
    times_per_mantis_pid_gff = []
    times_per_mantis_pid_gff_append = []

    with open(snakemake.output.gff, "w") as o_gff:
        for pid in gff.index:
            start_gff_pid = timer()
            # save (updated) annotations to output GFF
            pid_line = gff.loc[pid, :].copy()
            if pid in mantis.index:
                pid_line.loc["attr"] = ";".join(
                    [pid_line["attr"]] + list(mantis_links_generator(mantis.loc[pid, "Links"],
                                                                     sep=snakemake.params.sep)))
            o_gff.write("\t".join(pid_line.astype(str)) + "\n")
            end_gff_pid = timer()
            times_per_gff_pid.append(end_gff_pid - start_gff_pid)

            # save GFF attributes for featureCounts
            start_mantis_pid_total = timer()
            if pid in mantis.index:
                start_mantis_pid_links = timer()
                pid_links = mantis_links_pairs(mantis.loc[pid, "Links"], sep=snakemake.params.sep)
                end_mantis_pid_links = timer()
                times_per_mantis_pid_links.append(end_mantis_pid_links - start_mantis_pid_links)
                if pid_links is not None:
                    start_mantis_pid_gff_append = timer()
                    # save GFF info + new attribute to file
                    for key, key_df in pid_links.groupby("Key"):
                        with open(f"{snakemake.params.prefix}.{key}.gff", "a") as o_fc:
                            for key_v in key_df.Pair:
                                o_fc.write("\t".join(list(pid_line.drop(["attr"]).astype(str)) + [key_v]) + "\n")
                    end_mantis_pid_gff_append = timer()
                    times_per_mantis_pid_gff_append.append(end_mantis_pid_gff_append - start_mantis_pid_gff_append)

            end_mantis_pid_total = timer()
            times_per_mantis_pid_total.append(end_mantis_pid_total - start_mantis_pid_total)

            pids_processed += 1
            if pids_processed % 10000 == 0 or pids_processed == len(gff.index):
                sum_time_per_gff_pid = round(sum(times_per_gff_pid), 1)
                sum_time_per_mantis_pid_total = round(sum(times_per_mantis_pid_total), 1)
                sum_time_per_mantis_pid_links = round(sum(times_per_mantis_pid_links), 1)
                sum_time_per_mantis_pid_gff = round(sum(times_per_mantis_pid_gff), 1)
                sum_time_per_mantis_pid_gff_append = round(sum(times_per_mantis_pid_gff_append), 1)

                logging.info(f'{pids_processed} Mantis annotations processed. Sum times per: gff_pid: {sum_time_per_gff_pid},'
                             f' mantis_pid_total: {sum_time_per_mantis_pid_total}, mantis_pid_links: {sum_time_per_mantis_pid_links},'
                             f' mantis_pid_gff: {sum_time_per_mantis_pid_gff}, mantis_pid_gff_append: {sum_time_per_mantis_pid_gff_append}.')

                times_per_gff_pid = []
                times_per_mantis_pid_total = []
                times_per_mantis_pid_links = []
                times_per_mantis_pid_gff = []
                times_per_mantis_pid_gff_append = []

    logging.info('Done.')
