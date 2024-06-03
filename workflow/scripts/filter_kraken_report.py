###############
### IMPORTS ###
###############

import argparse

###############
### GLOBALS ###
###############

TAX_CODE_HIERARCHY = ('R', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S')

###################
### FUNCTIONS #####
###################

def parse_kraken_report_indent(report_file):
    taxa_dict = {}
    # Path stack will have tuples of (rank_code, tax_id, indent_level) as items
    path_stack = []

    indent_level = -1
    previous_indent_level = -1
    
    with open(report_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue

            level_info = {
                'percent_reads': float(parts[0]),
                'read_count': int(parts[1]),
                'read_count_direct': int(parts[2]),
                'minimizer_count': int(parts[3]),
                'unique_minimizer_count': int(parts[4]),
                'rank_code': parts[5],
                'tax_id': parts[6],
                'sci_name': parts[7]
            }

            # Get the indent level of the current line from level_info['sci_name']
            indent_level = len(level_info['sci_name']) - len(level_info['sci_name'].lstrip())
            level_info['indent_level'] = indent_level

            # Deal with special cases of levels unclassified (U) and root (R)
            if level_info['rank_code'] == 'U':
                # Add unclassified line at base outside of hierarchy
                taxa_dict['0'] = {'stats': level_info, 'children': None, 'parent': None}
                continue
            elif level_info['rank_code'] == 'R':
                taxa_dict['1'] = {'stats': level_info, 'children': [], 'parent': None}
                indent_level = 0
                path_stack.append(('R','1', indent_level))
                continue

            # Using the indent level, we can determine the hierarchy level
            if indent_level > previous_indent_level:
                # Going down the hierarchy, add item to parent taxon
                path_stack.append((level_info['rank_code'], level_info['tax_id'], indent_level))
            elif indent_level == previous_indent_level:
                # Same level, add item to parent taxon
                path_stack[-1] = (level_info['rank_code'], level_info['tax_id'], indent_level)
            else:
                # Go up stack until lower indent level is found, set level to parent hierarchy id
                while path_stack[-1][2] >= indent_level:
                    path_stack.pop()
                path_stack.append((level_info['rank_code'], level_info['tax_id'], indent_level))

            # Get the parent taxon id
            parent_tax_id = path_stack[-2][1]

            # Add taxon to dict
            taxa_dict[level_info['tax_id']] = {'stats': level_info, 'children': [], 'parent': parent_tax_id}

            # Add taxon to parent taxon
            taxa_dict[parent_tax_id]['children'].append(level_info['tax_id'])

            # Update previous level
            previous_indent_level = indent_level

            # print(f'name: {level_info['sci_name'].strip()}, rank_code: {level_info['rank_code']}, indent_level: {indent_level}, path_stack: {path_stack}')

    return taxa_dict


def get_all_indent_levels(taxa_dict):
    indent_levels = set()
    for key, value in taxa_dict.items():
        indent_levels.add(value['stats']['indent_level'])
    
    # Sort the indent levels
    indent_levels = sorted(list(indent_levels), reverse=True)

    return indent_levels


def filter_adjust_taxa(taxa_dict, min_reads=10, min_uniq_mini=1000):
    # Get all indent levels
    indent_levels = get_all_indent_levels(taxa_dict)

    # Track total reads before the ensure final read count is consistent
    total_read_prior = taxa_dict['1']['stats']['read_count'] + taxa_dict['0']['stats']['read_count']

    # Track the taxa to remove
    to_remove = []
    lineage_id = None

    for indent_level in indent_levels:
        for tax_id, tax_data in taxa_dict.items():
            tax_stats = tax_data['stats']
            if tax_stats['indent_level'] == indent_level:
                # We need to keep the root to always be able to calculate the total reads
                if tax_id not in ['0', '1'] and (tax_stats['read_count'] < min_reads or tax_stats['unique_minimizer_count'] < min_uniq_mini):
                    to_remove.append(tax_id)
                    parent_tax_id = tax_data['parent']
                    # Check if parent still exists in the dictionary
                    if parent_tax_id in taxa_dict:
                        taxa_dict[parent_tax_id]['children'].remove(tax_id)
                        # Now we need to go up until root and subtract the stats from entire lineage
                        lineage_id = parent_tax_id
                        while lineage_id:
                            lineage_tax_stats = taxa_dict[lineage_id]['stats']

                            lineage_tax_stats['read_count'] -= tax_stats['read_count']
                            # We dont want to subtract the direct read count, since its specific to the taxon and level
                            # parent_tax_stats['read_count_direct'] -= tax_stats['read_count_direct']
                            lineage_tax_stats['minimizer_count'] -= tax_stats['minimizer_count']
                            lineage_tax_stats['unique_minimizer_count'] -= tax_stats['unique_minimizer_count']

                            lineage_id = taxa_dict[lineage_id]['parent']
                        
                        # Add read count to the unclassified level
                        taxa_dict['0']['stats']['read_count'] += tax_stats['read_count']
                        taxa_dict['0']['stats']['read_count_direct'] += tax_stats['read_count']

    # Remove the taxa outside of the loop
    for tax_id in to_remove:
        if tax_id in taxa_dict:
            del taxa_dict[tax_id]

    # Finally recalculate the percentages
    total_reads = taxa_dict['1']['stats']['read_count'] + taxa_dict['0']['stats']['read_count']

    # Remove root if < min_reads
    if taxa_dict['1']['stats']['read_count'] < min_reads:
        del taxa_dict['1']

    # Check if the total reads are consistent
    assert total_reads == total_read_prior, f'Total reads before and after filtering do not match: {total_reads} != {total_read_prior}'

    for tax_id, tax_data in taxa_dict.items():
        tax_stats = tax_data['stats']
        tax_stats['percent_reads'] = (tax_stats['read_count'] / total_reads) * 100
    

def write_taxa(out_file_handle, taxa_dict, tax_id):
    # Initialize stack with tax_id
    stack = [tax_id]
    while stack:
        current_tax_id = stack.pop()
        tax_data = taxa_dict[current_tax_id]
        # Remove indent_level from stats dict
        del tax_data['stats']['indent_level']
        out_file_handle.write('\t'.join([str(val) if not isinstance(val, float) else f'{val:.6f}' for key, val in tax_data['stats'].items()]) + '\n')
        # Push children to the stack
        children_sorted = sorted(tax_data['children'], key=lambda x: taxa_dict[x]['stats']['sci_name'], reverse=True)
        for child in children_sorted:
            stack.append(child)

def write_filtered_report(taxa_dict, filename):
    with open(filename, 'w') as f:
        # Start with unclassified and root
        # Remove indent_level from stats dict
        del taxa_dict['0']['stats']['indent_level']
        f.write('\t'.join([str(val) if not isinstance(val, float) else f'{val:.6f}' for key, val in taxa_dict['0']['stats'].items()]) + '\n')

        # We only write if root exists (since if theres no root theres no taxa to write either)
        if '1' in taxa_dict:
            del taxa_dict['1']['stats']['indent_level']
            f.write('\t'.join([str(val) if not isinstance(val, float) else f'{val:.6f}' for key, val in taxa_dict['1']['stats'].items()]) + '\n')
        
            # Write the rest of the taxa with the correct hierarchy chaining from root's children and so on
            # Human sort taxa_dict['1']['children'] by taxa_dict['1']['stats']['sci_name']
            children_sorted = sorted(taxa_dict['1']['children'], key=lambda x: taxa_dict[x]['stats']['sci_name'], reverse=True)
            for child in children_sorted:
                write_taxa(f, taxa_dict, child)


# Example usage
# proj_dir = '/home/oss/Projects/Prospectomics/filt_kraken'
# kraken_report = f'{proj_dir}/mt.kraken.report'
# out_report = kraken_report.replace('.report', '_filtered.report')

# taxa_dict = parse_kraken_report_indent(kraken_report)

# filter_adjust_taxa(taxa_dict)
# write_filtered_report(taxa_dict, out_report)


if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description="Filter and adjust Kraken report file based on read count and unique minimizer count.")
    parser.add_argument("kraken_report", help="Path to the Kraken report file.")
    parser.add_argument("filtered_report", help="Path to the filtered Kraken report file.")
    parser.add_argument("--min_reads", help="Minimum read count to keep taxon.", default=10)
    parser.add_argument("--min_uniq_mini", help="Minimum unique minimizer count to keep taxon.", default=1000)
    args = parser.parse_args()

    min_reads = int(args.min_reads)
    min_uniq_mini = int(args.min_uniq_mini)

    taxa_dict = parse_kraken_report_indent(args.kraken_report)
    filter_adjust_taxa(taxa_dict, min_reads=min_reads, min_uniq_mini=min_uniq_mini)

    write_filtered_report(taxa_dict, args.filtered_report)