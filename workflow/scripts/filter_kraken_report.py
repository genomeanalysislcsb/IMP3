#!/usr/bin/env python3
"""
Filter and adjust Kraken 2 report files based on minimum read counts and unique minimizer thresholds.
"""

import argparse

# Hierarchical rank codes as defined by Kraken 2, though not strictly used in this script
TAX_CODE_HIERARCHY = ("R", "D", "K", "P", "C", "O", "F", "G", "S")


def parse_kraken_report_indent(report_file):
    """
    Parse the Kraken 2 report file preserving the indentation (hierarchical) structure.

    Parameters
    ----------
    report_file : str
        Path to the Kraken 2 report file.

    Returns
    -------
    dict
        Dictionary mapping tax_id (str) to a structure with:
        {
            "stats": {
                "percent_reads": float,
                "read_count": int,
                "read_count_direct": int,
                "minimizer_count": int,
                "unique_minimizer_count": int,
                "rank_code": str,
                "tax_id": str,
                "sci_name": str,
                "indent_level": int
            },
            "children": list of child tax_ids,
            "parent": parent tax_id or None
        }
    """
    taxa_dict = {}
    path_stack = []

    indent_level = -1
    previous_indent_level = -1

    with open(report_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if not parts or len(parts) < 8:
                continue

            level_info = {
                "percent_reads": float(parts[0]),
                "read_count": int(parts[1]),
                "read_count_direct": int(parts[2]),
                "minimizer_count": int(parts[3]),
                "unique_minimizer_count": int(parts[4]),
                "rank_code": parts[5],
                "tax_id": parts[6],
                "sci_name": parts[7],
            }

            # Determine indentation level from leading spaces in the scientific name
            indent_level = len(level_info["sci_name"]) - len(level_info["sci_name"].lstrip())
            level_info["indent_level"] = indent_level

            # Handle unclassified (U) and root (R)
            if level_info["rank_code"] == "U":
                taxa_dict["0"] = {
                    "stats": level_info,
                    "children": None,
                    "parent": None,
                }
                continue
            elif level_info["rank_code"] == "R":
                taxa_dict["1"] = {
                    "stats": level_info,
                    "children": [],
                    "parent": None,
                }
                indent_level = 0
                path_stack.append(("R", "1", indent_level))
                continue

            # Adjust the path stack based on indentation changes
            if indent_level > previous_indent_level:
                # Descending into a deeper level
                path_stack.append((level_info["rank_code"], level_info["tax_id"], indent_level))
            elif indent_level == previous_indent_level:
                # Same level as previous
                path_stack[-1] = (level_info["rank_code"], level_info["tax_id"], indent_level)
            else:
                # Ascending back up the hierarchy
                while path_stack and path_stack[-1][2] >= indent_level:
                    path_stack.pop()
                path_stack.append((level_info["rank_code"], level_info["tax_id"], indent_level))

            parent_tax_id = path_stack[-2][1] if len(path_stack) > 1 else None
            taxa_dict[level_info["tax_id"]] = {
                "stats": level_info,
                "children": [],
                "parent": parent_tax_id,
            }

            # Add the current tax_id to the parent's children list
            if parent_tax_id in taxa_dict and taxa_dict[parent_tax_id]["children"] is not None:
                taxa_dict[parent_tax_id]["children"].append(level_info["tax_id"])

            previous_indent_level = indent_level

    return taxa_dict


def filter_adjust_taxa(taxa_dict, min_reads=10, min_uniq_mini=1000):
    """
    Filters and adjusts the parsed Kraken 2 taxonomy dictionary based on thresholds.
    If root fails thresholds, all classified reads are moved to unclassified.
    
    Parameters
    ----------
    taxa_dict : dict
        Dictionary of taxa as produced by parse_kraken_report_indent().
    min_reads : int
        Minimum read count threshold.
    min_uniq_mini : int
        Minimum unique minimizer threshold.
    """
    # Identify root and unclassified keys
    unclassified_id = "0"
    root_id = "1"

    # Calculate total reads before filtering
    total_read_prior = 0
    if unclassified_id in taxa_dict:
        total_read_prior += taxa_dict[unclassified_id]["stats"]["read_count"]
    if root_id in taxa_dict:
        total_read_prior += taxa_dict[root_id]["stats"]["read_count"]

    # Check if root passes thresholds
    root_passes = False
    if root_id in taxa_dict:
        root_stats = taxa_dict[root_id]["stats"]
        root_passes = (root_stats["read_count"] >= min_reads and 
                       root_stats["unique_minimizer_count"] >= min_uniq_mini)
        
        # If root doesn't pass, transfer all its reads to unclassified and remove all taxa
        if not root_passes:
            # Ensure unclassified exists
            if unclassified_id not in taxa_dict:
                taxa_dict[unclassified_id] = {
                    "stats": {
                        "percent_reads": 0.0,
                        "read_count": 0,
                        "read_count_direct": 0,
                        "minimizer_count": 0,
                        "unique_minimizer_count": 0,
                        "rank_code": "U",
                        "tax_id": "0",
                        "sci_name": "unclassified",
                        "indent_level": 0
                    },
                    "children": None,
                    "parent": None,
                }
            
            # Transfer root reads to unclassified
            taxa_dict[unclassified_id]["stats"]["read_count"] += root_stats["read_count"]
            taxa_dict[unclassified_id]["stats"]["read_count_direct"] += root_stats["read_count_direct"]
            
            # Remove all taxa except unclassified
            taxa_keys = list(taxa_dict.keys())
            for key in taxa_keys:
                if key != unclassified_id:
                    del taxa_dict[key]
        else:
            # If root passes, filter individual taxa using existing logic
            # (Original filtering code would go here)
            to_remove = []
            queue = [root_id]
            while queue:
                current_id = queue.pop(0)
                if current_id not in [root_id, unclassified_id]:
                    stats = taxa_dict[current_id]["stats"]
                    if stats["read_count"] < min_reads or stats["unique_minimizer_count"] < min_uniq_mini:
                        to_remove.append(current_id)
                    else:
                        for c_id in taxa_dict[current_id]["children"]:
                            queue.append(c_id)
                else:
                    for c_id in taxa_dict[current_id]["children"]:
                        queue.append(c_id)
                        
            # Perform removal
            for node_id in to_remove:
                parent_id = taxa_dict[node_id]["parent"]
                if parent_id in taxa_dict and taxa_dict[parent_id]["children"]:
                    if node_id in taxa_dict[parent_id]["children"]:
                        taxa_dict[parent_id]["children"].remove(node_id)
                del taxa_dict[node_id]

    # Verify total reads remain consistent
    total_reads = 0
    if unclassified_id in taxa_dict:
        total_reads += taxa_dict[unclassified_id]["stats"]["read_count"]
    if root_id in taxa_dict:
        total_reads += taxa_dict[root_id]["stats"]["read_count"]

    assert (
        total_reads == total_read_prior
    ), f"Total reads before and after filtering do not match: {total_reads} != {total_read_prior}"

    # Recalculate percentages
    if total_reads > 0:
        for t_id, t_data in taxa_dict.items():
            t_data["stats"]["percent_reads"] = (t_data["stats"]["read_count"] / total_reads) * 100.0
    else:
        for t_id, t_data in taxa_dict.items():
            t_data["stats"]["percent_reads"] = 0.0


def write_taxa(out_file_handle, taxa_dict, tax_id):
    """
    Recursively writes a taxon and its descendants to an open file handle
    in a depth-first manner, sorting children by descending scientific name.

    Parameters
    ----------
    out_file_handle : file handle
        The file handle to write to.
    taxa_dict : dict
        Dictionary of taxa.
    tax_id : str
        Tax ID of the node to write (and its descendants).
    """
    stack = [tax_id]
    while stack:
        current_tax_id = stack.pop()
        tax_data = taxa_dict[current_tax_id]
        # Remove indent_level from output
        stats_copy = dict(tax_data["stats"])
        if "indent_level" in stats_copy:
            del stats_copy["indent_level"]

        # Convert stats into a tab-delimited line
        out_line_values = []
        for k in [
            "percent_reads",
            "read_count",
            "read_count_direct",
            "minimizer_count",
            "unique_minimizer_count",
            "rank_code",
            "tax_id",
            "sci_name",
        ]:
            val = stats_copy[k]
            if isinstance(val, float):
                out_line_values.append(f"{val:.6f}")
            else:
                out_line_values.append(str(val))
        out_file_handle.write("\t".join(out_line_values) + "\n")

        children_sorted = sorted(
            tax_data["children"],
            key=lambda x: taxa_dict[x]["stats"]["sci_name"],
            reverse=True,
        )
        for child in children_sorted:
            stack.append(child)


def write_filtered_report(taxa_dict, filename):
    """
    Writes the filtered Kraken 2 report to disk.

    Parameters
    ----------
    taxa_dict : dict
        Dictionary of taxa after filtering.
    filename : str
        Output file path.
    """
    with open(filename, "w") as f:
        # Write unclassified if present
        if "0" in taxa_dict:
            stats_copy = dict(taxa_dict["0"]["stats"])
            if "indent_level" in stats_copy:
                del stats_copy["indent_level"]
            out_line_values = []
            for k in [
                "percent_reads",
                "read_count",
                "read_count_direct",
                "minimizer_count",
                "unique_minimizer_count",
                "rank_code",
                "tax_id",
                "sci_name",
            ]:
                val = stats_copy[k]
                if isinstance(val, float):
                    out_line_values.append(f"{val:.6f}")
                else:
                    out_line_values.append(str(val))
            f.write("\t".join(out_line_values) + "\n")

        # Write root if present
        if "1" in taxa_dict:
            stats_copy = dict(taxa_dict["1"]["stats"])
            if "indent_level" in stats_copy:
                del stats_copy["indent_level"]
            out_line_values = []
            for k in [
                "percent_reads",
                "read_count",
                "read_count_direct",
                "minimizer_count",
                "unique_minimizer_count",
                "rank_code",
                "tax_id",
                "sci_name",
            ]:
                val = stats_copy[k]
                if isinstance(val, float):
                    out_line_values.append(f"{val:.6f}")
                else:
                    out_line_values.append(str(val))
            f.write("\t".join(out_line_values) + "\n")

            # Traverse children of root
            children_sorted = sorted(
                taxa_dict["1"]["children"],
                key=lambda x: taxa_dict[x]["stats"]["sci_name"],
                reverse=True,
            )
            for child in children_sorted:
                write_taxa(f, taxa_dict, child)


def main():
    """
    Command-line interface for filtering Kraken 2 report files
    based on read count and unique minimizer thresholds.
    """
    parser = argparse.ArgumentParser(
        description="Filter and adjust Kraken 2 report file based on read count "
                    "and unique minimizer thresholds. Implements a lenient "
                    "filtering paradigm."
    )
    parser.add_argument("kraken_report", help="Path to the Kraken report file.")
    parser.add_argument("filtered_report", help="Path to the filtered Kraken report file.")
    parser.add_argument(
        "--min_reads",
        type=int,
        default=10,
        help="Minimum read count to keep a taxon (default: 10).",
    )
    parser.add_argument(
        "--min_uniq_mini",
        type=int,
        default=1000,
        help="Minimum unique minimizer count to keep a taxon (default: 1000).",
    )
    args = parser.parse_args()

    # Parse the Kraken 2 report
    taxa_dict = parse_kraken_report_indent(args.kraken_report)

    # Filter the taxonomy
    filter_adjust_taxa(
        taxa_dict,
        min_reads=args.min_reads,
        min_uniq_mini=args.min_uniq_mini,
    )

    # Write out the filtered report
    write_filtered_report(taxa_dict, args.filtered_report)


if __name__ == "__main__":
    main()
