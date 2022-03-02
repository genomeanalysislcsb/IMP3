#!/usr/bin/env python3

##############################
# IMPORTS
##############################
import pandas

##############################
# FUNCTIONS
##############################
def read_in_gff(fname):
    """Read in a GFF"""
    # column names (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#description-of-the-format)
    gff_cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr"]
    # read in
    gff = pandas.read_csv(fname, sep="\t", header=None, names=gff_cols, comment="#")
    # set row names: feature ID extracted from the attributes column
    gff.set_index(keys=gff.attr.str.extract(r"(^|;)ID=(?P<id>\w+);")["id"], inplace=True, verify_integrity=True)
    return gff

def read_mantis_consensus(fname):
    """Read Mantis consensus table"""
    annot = pandas.read_csv(fname, sep="\t", header=0, index_col=0)
    return annot

def split_gff_attr(x, sep1=";", sep2="="):
    """Split keys/values found in a GFF attributes string"""
    if x and isinstance(x, str):
        return [y.split(sep2) for y in x.split(sep1)]
    return []

def split_mantis_links(x, sep1="|", sep2=":"):
    """Split keys/values found in a Mantis links string"""
    if x and isinstance(x, str):
        return [y.split(sep2) for y in x.split(sep1)]
    return []

def get_keys_from_gff_attr(x, sep1=";", sep2="="):
    """Get unique keys found in a GFF attributes string"""
    if x and isinstance(x, str):
        return set([y[0] for y in split_gff_attr(x=x, sep1=sep1, sep2=sep2)])
    return set()

def get_keys_from_mantis_links(x, sep1="|", sep2=":"):
    """Get unique keys found in a Mantis links string"""
    if x and isinstance(x, str):
        return set([y[0] for y in split_mantis_links(x=x, sep1=sep1, sep2=sep2)])
    return set()

def mantis_links2dataframe(links_str, sep1="|", sep2=":"):
    """
    Create a dataframe of keys and values from Mantis consensus annotation (links) string

    Parameters
    ----------
    links_str : str
        Links string provided by Mantis in the consensus output table
    sep : str
        Character/string used to separate anotations

    Return
    ------
    pandas.DataFrame
        Table: Key, Value
    """
    # dataframe of keys and values
    df_links = pandas.DataFrame(split_mantis_links(x=links_str, sep1=sep1, sep2=sep2), columns=["Key", "Value"])
    if df_links.shape[0] == 0:
        return None
    # remove rows: if key is "description" (irrelevant for GFF)
    df_links = df_links.loc[df_links.Key != "description",:]
    # add column: key/value pairs (<key>=<value>)
    df_links = df_links.assign(Pair = df_links.apply(func="=".join, axis=1))
    return df_links

def mantis_gff_attrs_generator(links_str, sep1="|", sep2=":"):
    """
    Generator yielding summarized Mantis annotations in GFF attribute format

    Parameters
    ----------
    links_str : str
        Links string provided by Mantis in the consensus output table
    sep1 : str
        Character/string used to separate anotations
    sep2 : str
        Character/string used to separate key and value in anotations

    Yields
    ------
    str
        Summarized annotation for a particular key in GFF attribute format, i.e. <key>=<value1,value2,...>
    """
    # data frame of keys and values
    df_links = mantis_links2dataframe(links_str, sep1=sep1, sep2=sep2)
    # generator
    if df_links is not None:
        for gr_name, gr_df in df_links.groupby("Key"):
            yield "{}={}".format(gr_name, ",".join(gr_df.Value.astype(str)))

##############################
# MAIN
##############################

if __name__ == "__main__":
    # read in
    gff    = read_in_gff(snakemake.input.gff)
    mantis = read_mantis_consensus(snakemake.input.mantis)

    # check attribute keys appearing in given input files: there should be no overlap
    gff_keys = gff.attr.apply(func=get_keys_from_gff_attr).values
    gff_keys = set.union(*gff_keys)
    mantis_keys = mantis.Links.apply(func=get_keys_from_mantis_links, sep1=snakemake.params.sep).values
    mantis_keys = set.union(*mantis_keys)
    assert len(gff_keys.intersection(mantis_keys)) == 0, "There are common attribute keys in {}: {}".format(input, gff_keys.intersection(mantis_keys))

    opend_files = set()
    for pid in mantis.index:
        # add GFF attributes based on Mantis annotations
        gff.loc[pid, "attr"] = ";".join(
            [gff.loc[pid, "attr"]] + list(mantis_gff_attrs_generator(mantis.loc[pid, "Links"], sep1=snakemake.params.sep))
        )

        # add gff attributes for featureCounts
        curlinks = mantis_links2dataframe(mantis.loc[pid, "Links"], sep1=snakemake.params.sep)
        if curlinks is not None:
            curlinks.drop(columns="Value", inplace=True)
            curline = gff.loc[[pid]*len(curlinks.index),gff.columns != "attr"].copy()
            curline = curline.set_index(curlinks.index[:len(curline)]).reindex(curlinks.index, method="ffill")
            curall = pandas.concat([curline,curlinks], axis=1)
            for k, k_links in curall.groupby("Key"):
                op = snakemake.params.prefix + "." + k + ".gff"
                k_links.drop(columns="Key").to_csv(op, sep="\t", header=False, index=False, mode="a" if op in opend_files else "w")
                opend_files.add(op)
    
    # save
    gff.to_csv(snakemake.output.gff, sep="\t", header=False, index=False)
