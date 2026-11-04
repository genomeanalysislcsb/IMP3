#!/usr/bin/env python3

##############################
# IMPORTS
##############################
import pandas
import logging
import os

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

    if not df_links.empty:
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
    logging.info('Reading annotation gff.')
    gff = read_in_gff(snakemake.input.gff)
    logging.info('Reading Mantis output.')
    mantis = read_mantis_consensus(snakemake.input.mantis)

    # # Fill NAs in Links col to prevent breaking of processing
    # (assignment, not inplace=True on the column: with pandas' copy-on-write, the default from pandas 3, an
    # inplace fillna on mantis['Links'] does not change mantis)
    mantis['Links'] = mantis['Links'].fillna('Annotation_status:Unannotated')
    # print(mantis[mantis['Links'] == 'Unannotated'])

    # check attribute keys appearing in given input files: there should be no overlap
    logging.info('Assure no overlap in input file keys.')
    gff_keys = gff.attr.apply(lambda x: set([y.split("=")[0] for y in x.split(";")])).values
    gff_keys = set.union(*gff_keys)
    # mantis_keys = mantis.Links.apply(lambda x: set([y.split(":")[0] for y in x.split(snakemake.params.sep)])).values
    mantis_keys = mantis.Links.apply(lambda x: set([y.split(":")[0] for y in x.split(snakemake.params.sep)])
                                               if isinstance(x, str) else set()).values

    mantis_keys = set.union(*mantis_keys)
    assert len(gff_keys.intersection(mantis_keys)) == 0,\
           f'There are common attribute keys in {input}: {gff_keys.intersection(mantis_keys)}'

    # add GFF attributes based on Mantis' annotations
    logging.info('Add GFF attributes based on Mantis annotations and GFF attributes per db for featureCounts.')

    # Print keys to write
    logging.info(f'Will write gffs for refs: {mantis_keys}')

    pids_processed = 0

    # Ensure per ref gff dir is made
    if not os.path.exists(os.path.dirname(snakemake.params.prefix)):
        os.makedirs(os.path.dirname(snakemake.params.prefix))

    with open(snakemake.output.gff, "w") as o_gff:
        n_entries = len(gff.index)
        for pid in gff.index:
            # Save (updated) annotations to output GFF
            pid_line = gff.loc[pid, :].copy()
            if pid in mantis.index:
                pid_line.loc["attr"] = ";".join(
                    [pid_line["attr"]] + list(mantis_links_generator(mantis.loc[pid, "Links"],
                                                                     sep=snakemake.params.sep)))
            o_gff.write("\t".join(pid_line.astype(str)) + "\n")

            # Save GFF attributes for featureCounts
            if pid in mantis.index:
                pid_links = mantis_links_pairs(mantis.loc[pid, "Links"], sep=snakemake.params.sep)
                if pid_links is not None:
                    # Save GFF info + new attribute to file
                    for key, key_df in pid_links.groupby("Key"):
                        with open(f"{snakemake.params.prefix}.{key}.gff", "a") as o_fc:
                            for key_v in key_df.Pair:
                                o_fc.write("\t".join(list(pid_line.drop(["attr"]).astype(str)) + [key_v]) + "\n")

            pids_processed += 1
            percent_processed = round(pids_processed / n_entries * 100, 1)
            
            if pids_processed % 10000 == 0 or pids_processed == n_entries:
                logging.info(f'{pids_processed}/{n_entries} annotations processed ({percent_processed}%)')

    logging.info('Done.')
