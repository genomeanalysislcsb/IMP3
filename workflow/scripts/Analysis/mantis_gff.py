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
    # read in
    annot = pandas.read_csv(fname, sep="\t", header=0, index_col=0)
    return annot

def mantis_gff_attrs_generator(links_str, sep="|"):
    """
    Generator yielding summarized Mantis' annotations in GFF attribute format

    Parameters
    ----------
    links_str : str
        Links string provided by Mantis in the consensus output table
    sep : str
        Character/string used to separate anotations

    Yields
    ------
    str
        Summarized annotation for a particular key in GFF attribute format, i.e. <key>=<value1,value2,...>
    """
    # data frame of keys and values
    df_links = pandas.DataFrame([tuple(hit.split(":")) for hit in links_str.split(sep)], columns=["Key","Value"])
    # remove description
    df_links = df_links.loc[df_links.Key != "description",:]
    # generator
    for gr_name, gr_df in df_links.groupby("Key"):
        yield "{}={}".format(gr_name, ",".join(gr_df.Value.astype(str)))

##############################
# MAIN
##############################

if __name__ == '__main__':
    # read in
    gff    = read_in_gff(snakemake.input.gff)
    mantis = read_mantis_consensus(snakemake.input.mantis)

    # check attribute keys appearing in given input files: there should be no overlap
    gff_keys = gff.attr.apply(lambda x: set([y.split("=")[0] for y in x.split(";")])).values
    gff_keys = set.union(*gff_keys)
    mantis_keys = mantis.Links.apply(lambda x: set([y.split(":")[0] for y in x.split(snakemake.params.sep)])).values
    mantis_keys = set.union(*mantis_keys)
    assert len(gff_keys.intersection(mantis_keys)) == 0, "There are common attribute keys in {}: {}".format(input, gff_keys.intersection(mantis_keys))

    # add GFF attributes based on Mantis' annotations
    for pid in mantis.index:
        gff.loc[pid, "attr"] = ";".join(
            [gff.loc[pid, "attr"]] + list(mantis_gff_attrs_generator(mantis.loc[pid, "Links"], sep=snakemake.params.sep))
        )
    
    # save
    gff.to_csv(snakemake.output.gff, sep="\t", header=False, index=False)