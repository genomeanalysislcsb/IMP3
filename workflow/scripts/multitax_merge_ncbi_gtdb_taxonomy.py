#!/usr/bin/env python3
#
# multitax_merge_ncbi_gtdb_taxonomy.py
#
# This script merges NCBI and GTDB taxonomies into a unified taxdump format
# (nodes.dmp, names.dmp, merged.dmp). It prunes the prokaryotic branches from
# NCBI and grafts the GTDB taxonomy in their place.
#
# It is designed to be called from a Snakemake rule and accepts all
# necessary paths as command-line arguments.

import sys
import os
import gzip
import re
import argparse
from collections import defaultdict

# Import necessary multitax classes
try:
    from multitax import NcbiTx, GtdbTx
except ImportError:
    print("Error: multitax library not found. Please install it (`pip install multitax`)", file=sys.stderr)
    sys.exit(1)

# --- Configuration ---
# Output filenames (will be placed in the specified output directory)
MERGED_NODES_FILENAME = "nodes.dmp"
MERGED_NAMES_FILENAME = "names.dmp"
MERGED_MERGED_FILENAME = "merged.dmp" # We need to write a merged.dmp for completeness

# --- Taxonomy Settings ---
NCBI_BACTERIA_ID = "2"
NCBI_ARCHAEA_ID = "2157"
GTDB_START_ID = 3000000000 # Starting ID for new numeric IDs for GTDB nodes

# GTDB Rank prefixes and full names (order matters!)
GTDB_RANK_MAP = {
    "d__": "domain", "p__": "phylum", "c__": "class", "o__": "order",
    "f__": "family", "g__": "genus", "s__": "species"
}
GTDB_RANK_NAMES_ORDERED = list(GTDB_RANK_MAP.values())
VALID_RANK_NAMES = set(GTDB_RANK_NAMES_ORDERED)
RANK_ORDER_MAP = {rank: i for i, rank in enumerate(GTDB_RANK_NAMES_ORDERED)}
UNKNOWN_RANK_ORDER = len(RANK_ORDER_MAP)


# --- Helper Functions ---
def print_stderr(message, **kwargs):
    """Prints messages to stderr."""
    print(message, file=sys.stderr)

def print_debug(message, debug_flag):
    """Prints debug messages to stderr only if debug_flag is True."""
    if debug_flag:
        print(f"DEBUG: {message}", file=sys.stderr)

def parse_gtdb_node_string(node_str):
    """ Parses 'p__Proteobacteria' into ('p__', 'Proteobacteria', 'phylum'). """
    for prefix, rank_name in GTDB_RANK_MAP.items():
        if node_str.startswith(prefix):
            name_part = node_str[len(prefix):]
            if not name_part:
                return prefix, None, rank_name
            return prefix, name_part, rank_name
    return None, node_str, None

def get_corrected_gtdb_name(node_id, node_name, node_rank, lineage_names, debug_flag=False):
    """
    Corrects GTDB names that are 'Unknown' or 'Bin<number>' based on the
    last known valid parent rank name from the provided lineage.
    """
    if node_name is None or node_name == "": return node_name

    original_name = node_name
    is_unknown = original_name.lower() == "unknown"
    is_bin = bool(re.fullmatch(r"Bin\d+", original_name, re.IGNORECASE))

    if not is_unknown and not is_bin: return original_name
    if not node_rank or node_rank == "no rank":
        print_debug(f"Correction Skipped: Cannot correct '{original_name}' node '{node_id}' due to missing rank.", debug_flag)
        return original_name

    last_valid_name = next((name for name in reversed(lineage_names) if name and name.strip()), None)

    if last_valid_name is None:
        print_debug(f"Correction Warning: No valid ancestor for '{node_id}'. Using original name.", debug_flag)
        return original_name

    safe_node_rank = node_rank.replace(" ", "_").lower()
    base_parent_name = last_valid_name.split('_')[0] # Heuristic to get a cleaner base name

    if is_unknown:
        corrected_name = f"Unknown_{base_parent_name}_{safe_node_rank}"
        print_debug(f"Correction (Unknown): Node='{node_id}', Parent='{last_valid_name}' -> Corrected='{corrected_name}'", debug_flag)
        return corrected_name
    elif is_bin:
        corrected_name = f"{base_parent_name}_{safe_node_rank}_{original_name}"
        print_debug(f"Correction (Bin): Node='{node_id}', Parent='{last_valid_name}' -> Corrected='{corrected_name}'", debug_flag)
        return corrected_name
    return original_name

def get_rank_sort_key(node_id, combined_nodes_info):
    """ Returns a sort key tuple (rank_order, node_id) for sorting nodes top-down. """
    rank = combined_nodes_info.get(node_id, {}).get('rank', 'no rank')
    rank_order = RANK_ORDER_MAP.get(rank.lower() if rank else None, UNKNOWN_RANK_ORDER)
    return (rank_order, node_id)


# --- Main Function ---
def main(args):
    """ Main execution function """
    DEBUG = args.debug

    # 1. Validate paths and create output directory
    print_stderr("Step 1/5: Initializing and validating paths...")
    os.makedirs(args.output_dir, exist_ok=True)
    
    ncbi_files = [
        os.path.join(args.ncbi_taxdump_dir, "nodes.dmp"),
        os.path.join(args.ncbi_taxdump_dir, "names.dmp"),
        os.path.join(args.ncbi_taxdump_dir, "merged.dmp")
    ]
    for f in ncbi_files + args.gtdb_files:
        if not os.path.exists(f):
            print_stderr(f"Error: Required input file not found: {f}")
            sys.exit(1)

    # 2. Load and prepare NCBI taxonomy
    print_stderr("Step 2/5: Loading and pruning NCBI taxonomy...")
    try:
        ncbi_tax = NcbiTx(files=ncbi_files, root_node="1")
        print_stderr(f" -> NCBI loaded with {len(ncbi_tax._nodes)} nodes.")
        
        # Prune the original prokaryotic branches
        nodes_to_prune = [NCBI_BACTERIA_ID, NCBI_ARCHAEA_ID]
        ncbi_tax.prune([n for n in nodes_to_prune if n in ncbi_tax._nodes])
        print_stderr(" -> Pruning of NCBI prokaryotes complete.")

        # Re-add the domain nodes, linked to root
        ncbi_tax._nodes[NCBI_BACTERIA_ID] = ncbi_tax.root_node
        ncbi_tax._ranks[NCBI_BACTERIA_ID] = "domain"
        ncbi_tax._names[NCBI_BACTERIA_ID] = "Bacteria"
        
        ncbi_tax._nodes[NCBI_ARCHAEA_ID] = ncbi_tax.root_node
        ncbi_tax._ranks[NCBI_ARCHAEA_ID] = "domain"
        ncbi_tax._names[NCBI_ARCHAEA_ID] = "Archaea"
        
        ncbi_tax._reset_aux_data()
        print_stderr(" -> NCBI Bacteria and Archaea domains re-linked to root.")

    except Exception as e:
        print_stderr(f"Error during NCBI taxonomy processing: {e}")
        sys.exit(1)

    # 3. Process GTDB taxonomy to build an intermediate structure
    print_stderr("Step 3/5: Processing GTDB taxonomy files...")
    combined_nodes_info = {}
    gtdb_root = '1' # GTDB lineages are relative to the root

    for tax_file in args.gtdb_files:
        print_stderr(f" -> Reading {os.path.basename(tax_file)}")
        with gzip.open(tax_file, 'rt') as f:
            for line in f:
                try:
                    _, lineage_str = line.strip().split('\t')
                    lineage_parts = [part.strip() for part in lineage_str.split(';') if part.strip()]
                    if not lineage_parts: continue

                    parent_lineage_key = gtdb_root
                    for i, node_str in enumerate(lineage_parts):
                        current_lineage_key = ";".join(lineage_parts[:i+1])
                        if current_lineage_key in combined_nodes_info:
                            parent_lineage_key = current_lineage_key
                            continue
                        
                        prefix, name_part, rank_name = parse_gtdb_node_string(node_str)
                        if not rank_name: continue

                        combined_nodes_info[current_lineage_key] = {
                            'parent_lineage_key': parent_lineage_key,
                            'rank': rank_name,
                            'original_name': name_part,
                            'gtdb_node_id': node_str
                        }
                        parent_lineage_key = current_lineage_key
                except ValueError:
                    continue # Skip malformed lines

    # 4. Correct names, assign numeric IDs, and add to the main taxonomy object
    print_stderr(f"Step 4/5: Integrating {len(combined_nodes_info)} GTDB nodes...")
    nodes_to_process_sorted = sorted(combined_nodes_info.keys(), key=lambda k: get_rank_sort_key(k, combined_nodes_info))
    
    corrected_names_map = {}
    lineage_key_to_numeric_id = {}
    next_gtdb_numeric_id = GTDB_START_ID

    for lineage_key in nodes_to_process_sorted:
        node_info = combined_nodes_info[lineage_key]
        
        # Correct name
        parent_lineage_key = node_info['parent_lineage_key']
        ancestor_keys = lineage_key.split(';')[:-1]
        ancestor_names = [corrected_names_map.get(k) for k in ancestor_keys]
        corrected_name = get_corrected_gtdb_name(node_info['gtdb_node_id'], node_info['original_name'], node_info['rank'], ancestor_names, DEBUG)
        corrected_names_map[lineage_key] = corrected_name

        # Assign numeric ID
        numeric_node_id = str(next_gtdb_numeric_id)
        if lineage_key == "d__Bacteria": numeric_node_id = NCBI_BACTERIA_ID
        elif lineage_key == "d__Archaea": numeric_node_id = NCBI_ARCHAEA_ID
        else: next_gtdb_numeric_id += 1
        lineage_key_to_numeric_id[lineage_key] = numeric_node_id

        # Add to taxonomy object (if not a domain)
        if numeric_node_id not in (NCBI_BACTERIA_ID, NCBI_ARCHAEA_ID):
            parent_numeric_id = lineage_key_to_numeric_id.get(parent_lineage_key)
            if not parent_numeric_id:
                # Determine domain parent for top-level phyla
                if lineage_key.startswith("d__Bacteria;"): parent_numeric_id = NCBI_BACTERIA_ID
                elif lineage_key.startswith("d__Archaea;"): parent_numeric_id = NCBI_ARCHAEA_ID
                else: parent_numeric_id = ncbi_tax.root_node # Fallback
            
            if corrected_name and parent_numeric_id:
                ncbi_tax._nodes[numeric_node_id] = parent_numeric_id
                ncbi_tax._ranks[numeric_node_id] = node_info['rank']
                ncbi_tax._names[numeric_node_id] = corrected_name

    ncbi_tax._reset_aux_data()
    print_stderr(f" -> Integration complete. Final taxonomy has {len(ncbi_tax._nodes)} nodes.")

    # 5. Write merged taxonomy to DMP files
    print_stderr(f"Step 5/5: Writing final taxdump files to {args.output_dir}...")
    try:
        nodes_out_path = os.path.join(args.output_dir, MERGED_NODES_FILENAME)
        names_out_path = os.path.join(args.output_dir, MERGED_NAMES_FILENAME)
        merged_out_path = os.path.join(args.output_dir, MERGED_MERGED_FILENAME)

        with open(nodes_out_path, "w") as fnodes, open(names_out_path, "w") as fnames:
            sorted_node_ids = sorted(ncbi_tax._nodes.keys(), key=lambda x: int(x))
            for node_id in sorted_node_ids:
                node_str = str(node_id)
                parent_str = str(ncbi_tax.parent(node_str))
                rank_str = ncbi_tax.rank(node_str)
                name_str = ncbi_tax.name(node_str).replace('\t', ' ').replace('\n', ' ')

                fnodes.write(f"{node_str}\t|\t{parent_str}\t|\t{rank_str}\t|\t\t|\n")
                fnames.write(f"{node_str}\t|\t{name_str}\t|\t\t|\tscientific name\t|\n")
        
        # Write the original NCBI merged.dmp content to the new location
        with open(ncbi_files[2], 'r') as f_in, open(merged_out_path, 'w') as f_out:
            f_out.write(f_in.read())

        print_stderr(" -> Taxdump files written successfully.")

    except Exception as e:
        print_stderr(f"Error writing final DMP files: {e}")
        sys.exit(1)

    print_stderr("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge NCBI and GTDB taxonomies into a unified taxdump.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--ncbi-taxdump-dir",
        required=True,
        help="Path to the directory containing NCBI nodes.dmp, names.dmp, and merged.dmp."
    )
    parser.add_argument(
        "--gtdb-files",
        nargs='+',
        required=True,
        help="One or more paths to GTDB taxonomy files (e.g., ar_taxonomy.tsv.gz bac_taxonomy.tsv.gz)."
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where the merged taxdump files will be written."
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable detailed debug printing to stderr."
    )
    args = parser.parse_args()
    main(args)
