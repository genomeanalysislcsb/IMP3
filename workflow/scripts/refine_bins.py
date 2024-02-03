#!/usr/bin/env python3

import os
import itertools
import argparse
import subprocess
from multiprocessing import Pool
from Bio import SeqIO
# import hdbscan
# from clustpy.partition import UniDip
from scipy.spatial.distance import squareform
import umap
import matplotlib.pyplot as plt
import seaborn as sns


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

def parse_arguments():
    parser = argparse.ArgumentParser(description="Refine bin sets using binning_refiner, assess quality using CheckM2 and consolidate them into a final set of bins.")
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help="output directory")
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help="number of threads (default=1)")
    # parser.add_argument('-m', '--memory', type=int, default=40,
    #                     help="memory available in GB (default=40)")
    parser.add_argument('-c', '--min_completion', type=int, default=70,
                        help="minimum completion of bins (default=70)")
    parser.add_argument('-x', '--max_contamination', type=int, default=10,
                        help="maximum contamination of bins that is acceptable (default=10)")
    parser.add_argument('-i', '--bin_folders', action=SplitArgs, required=True,
                        help="comma-separated list of bin folders to use for refinement")
    parser.add_argument('--skip-refinement', action='store_true', default=False,
                        help="don't use binning_refiner to come up with refined bins based on combinations of binner outputs")
    parser.add_argument('--skip-checkm', action='store_true', default=False,
                        help="don't run CheckM to assess bins")
    parser.add_argument('--skip-consolidation', action='store_true', default=True,
                        help="choose the best version of each bin from all bin refinement iteration")
    parser.add_argument('--keep-ambiguous', action='store_true', default=False,
                        help="for contigs that end up in more than one bin, keep them in all bins (default: keeps them only in the best bin)")
    parser.add_argument('--remove-ambiguous', action='store_true', default=False,
                        help="for contigs that end up in more than one bin, remove them in all bins (default: keeps them only in the best bin)")
    parser.add_argument('--quick', action='store_true', default=True,
                        help="adds --reduced_tree option to checkm, reducing runtime, especially with low memory")
    # get binning refiner path, default is the same dir as this script
    parser.add_argument('--binning_refiner_path', type=str, default=os.path.dirname(os.path.realpath(__file__)),
                        help="path to binning_refiner script")
    parser.add_argument('--checkm_db', type=str,
                        help="path to checkm database")

    args = parser.parse_args()
    return args


def create_symlinks_and_prepare_combinations(bin_folders, output_dir, include_singles=False):
    n = len(bin_folders)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    combinations = []
    # Generate combinations for all sizes from n-1 down to 2
    for size in range(n-1, 1, -1):  # Start from n-1, end at 2 (inclusive)
        for i, combination in enumerate(itertools.combinations(bin_folders, size), start=1):
            # Combination name based on input bin folder names
            binner_names = [os.path.basename(binner_folder) for binner_folder in combination]
            subdir_name = f"combination_size{size}_{i}_{'_'.join(binner_names)}"
            subdir_path = os.path.join(output_dir, subdir_name)
            os.makedirs(subdir_path, exist_ok=True)
            
            for bin_folder in combination:
                symlink_target = os.path.abspath(bin_folder)
                symlink_name = os.path.join(subdir_path, os.path.basename(bin_folder))
                if not os.path.exists(symlink_name):
                    os.symlink(symlink_target, symlink_name)
            
            combinations.append({'index': f'combination_size{size}_{i}', 'dir': subdir_path})

    # Optionally include original inputs as single entities
    if include_singles:
        for j, bin_folder in enumerate(bin_folders, start=1):
            binner_name = os.path.basename(bin_folder)
            subdir_name = f"single_{j}_{binner_name}"
            subdir_path = os.path.join(output_dir, subdir_name)
            os.makedirs(subdir_path, exist_ok=True)
            
            symlink_target = os.path.abspath(bin_folder)
            symlink_name = os.path.join(subdir_path, os.path.basename(bin_folder))
            if not os.path.exists(symlink_name):
                os.symlink(symlink_target, symlink_name)
            
            combinations.append({'index': f'single_{j}', 'dir': subdir_path})

    return combinations



def run_binning_refiner(combination, output_dir, refiner_path=''):
    # Construct binning_refiner command
    output_prefix = os.path.join(output_dir, f"{combination['index']}_refined")
    print(f"Running binning_refiner on {combination['dir']}")
    # raise Exception
    command = [
        f"{refiner_path}/binning_refiner.py",
        "-i", combination['dir'],
        "-o", output_prefix,
        "-p", output_prefix,
        "-m", "512",
        "-q", "False"
    ]

    # Run binning_refiner
    subprocess.run(command, check=True)

    refined_out = os.path.join(output_prefix, 'Refined')

    # Rename each bin fasta in the output directory to include the combination index
    for bin_file in [f for f in os.listdir(refined_out) if f.endswith(('.fa', '.fasta', '.fna'))]:
        # Strip extension and replace with `.fasta`
        bin_file = os.path.splitext(bin_file)[0] + '.fasta'
        os.rename(os.path.join(refined_out, bin_file), os.path.join(refined_out, f"{combination['index']}_{bin_file}"))

    # Return output prefix to be used as input for CheckM2
    return refined_out


def run_checkm2(input_folder, output_folder, db_path, threads):
    # Construct CheckM2 command
    command = [
        "checkm2",
        "predict",
        "--threads", str(threads),
        "--input", input_folder,
        "--output-directory", output_folder,
        "--extension", "fasta",
        "--force",
        "--database_path", db_path,
    ]

    # Run CheckM2
    subprocess.run(command, check=True)
    print(f"Completed checkm2 prediction for {input_folder}")


def run_checkm2_on_refined_outputs(refined_outputs, checkm_db, threads):
    binning_checkm_stats = {}
    for refined_output in refined_outputs:
        input_folder = f"{refined_output}"
        output_folder = f"{input_folder}_checkm2"
        # Ensure the output directory exists
        os.makedirs(output_folder, exist_ok=True)
        run_checkm2(input_folder, output_folder, checkm_db, threads)

        # Parse CheckM2 output
        checkm2_output_file = f"{output_folder}/quality_report.tsv"

        ref_out_name = os.path.basename(refined_output)
        binning_checkm_stats[ref_out_name] = parse_checkm2_output(checkm2_output_file)
    return binning_checkm_stats


def get_bin_cont_lengths(binner_out_paths):
    bin_cont_lenghts = {}
    for binner_out_path in binner_out_paths:
        # Get binning output name
        bin_set = os.path.basename(binner_out_path)
        bin_cont_lenghts[bin_set] = {}
        # Get bin files
        bin_files = [f for f in os.listdir(binner_out_path) if f.endswith(('.fa', '.fasta', '.fna'))]
        for bin_file in bin_files:
            bin_cont_lenghts[bin_set][bin_file] = {}
            # Full bin path
            bin_file_path = os.path.join(binner_out_path, bin_file)
            for contig in SeqIO.parse(bin_file_path, 'fasta'):
                cont_length = len(contig)
                bin_cont_lenghts[bin_set][bin_file][contig.id] = cont_length
    return bin_cont_lenghts


def parse_checkm2_output(checkm2_output_file, return_header=False):
    # Fetch info from `quality_report.tsv`
    metrics_dict = {}
    with open(checkm2_output_file, 'r') as f:
        # Fetch header
        for line in f:
            if line.startswith('Bin Id'):
                # Ignore first column
                header = line.strip().split('\t')[1:]
                break
        # Fetch bin metrics
        for line in f:
            if line.startswith('Bin Id'):
                continue
            cut = line.strip().split('\t')
            bin_id = cut[0]
            metrics_dict[bin_id] = {}
            for i, header_item in enumerate(header):
                metrics_dict[bin_id][header_item] = cut[i]
    if return_header:
        header = ['Bin Id'] + header
        return metrics_dict, header
    else:
        return metrics_dict


def calculate_weighted_bin_similarity(bin_cont_lenghts):
    similarity_scores = {}
    
    # Calculate similarity between each pair of bins from each set
    for bin_set_i_name, bin_set_i in bin_cont_lenghts.items():
        for bin_set_j_name, bin_set_j in bin_cont_lenghts.items():
            # Dont compare a bin set to itself
            if bin_set_i_name != bin_set_j_name:
                for bin_i, contigs_i in bin_set_i.items():
                    for bin_j, contigs_j in bin_set_j.items():
                        # Get bin sizes and contig names
                        bin_i_size = sum(contigs_i.values())
                        bin_i_contigs = set(contigs_i.keys())
                        bin_j_size = sum(contigs_j.values())
                        bin_j_contigs = set(contigs_j.keys())

                        # Find shared contigs and calculate shared length
                        shared_length = 0
                        shared_contigs = bin_i_contigs.intersection(bin_j_contigs)
                        shared_length += sum(contigs_i[contig] for contig in shared_contigs)

                        # Calculate similarity
                        similarity = 2 * shared_length / (bin_i_size + bin_j_size)
                        similarity_scores[(f"{bin_set_i_name}:{bin_i}", f"{bin_set_j_name}:{bin_j}")] = similarity
                
    return similarity_scores


def plot_clusters(distance_matrix, clusters, filename='umap_clusters.svg'):
    # Initialize UMAP
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='precomputed')
    
    # Fit and transform the distance matrix to 2D
    embedding = reducer.fit_transform(distance_matrix)
    
    # Plot the embedding with cluster labels
    plt.figure(figsize=(12, 10))
    sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=clusters, palette='Spectral', s=50, alpha=0.8, legend='full')
    plt.title('Clusters visualized with UMAP')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.legend(title='Cluster')

    # Save the plot to a file
    plt.savefig(filename, format='svg', dpi=300)  # Adjust 'format' as needed (e.g., 'pdf', 'eps')
    plt.close()  # Close the plot to free memory


def choose_best_bin_per_cluster(bin2cluster, checkm_stats, min_completion=60, max_contamination=10):
    # Initialize a dict to store the best bin in each cluster
    best_bin_per_cluster = {}
    cluster_scores = {}

    for bin_id, cluster_id in bin2cluster.items():
        bin_set, bin_name = bin_id.split(':')
        bin_stats = checkm_stats[bin_set][bin_name]

        # Dont consider bins that dont meet the minimum completion and maximum contamination criteria
        if bin_stats['completeness'] < min_completion or bin_stats['contamination'] > max_contamination:
            continue

        # Calculate score
        score = bin_stats['completeness'] - 2 * bin_stats['contamination']

        # Update best bin if this bin has a higher score than the current best, or if this cluster has no best bin yet
        if cluster_id not in cluster_scores or score > cluster_scores[cluster_id]:
            cluster_scores[cluster_id] = score
            # Keep bin with its stats
            best_bin_per_cluster[cluster_id] =  bin_id # (bin_id, bin_stats)
        # If the score is the same, print a warning
        elif score == cluster_scores[cluster_id]:
            print(f"Warning: bin {bin_id} has the same score as {best_bin_per_cluster[cluster_id]}")

    return best_bin_per_cluster


def identify_and_resolve_ambiguous_contigs(bin_sets, best_bin_per_cluster, keep_ambiguous, remove_ambiguous):
    contig_to_bins = {}
    # Aggregate contig appearances
    for cluster_id, best_bin in best_bin_per_cluster.items():
        # bin_set_name, checkm_stats = best_bin
        bin_set, bin_name = best_bin.split(':')
        bin_contigs = bin_sets[bin_set][bin_name]
        for contig in bin_contigs:
            if contig not in contig_to_bins:
                contig_to_bins[contig] = [best_bin]
            else:
                contig_to_bins[contig].append(best_bin)

    # Identify ambiguous contigs
    ambiguous_contigs = {contig: bins for contig, bins in contig_to_bins.items() if len(bins) > 1}

    # log ambiguous contigs
    print(f"Found {len(ambiguous_contigs)} ambiguous contigs")

    # Resolve ambiguity
    if remove_ambiguous:
        # Example resolution strategy: remove ambiguous contigs from all but the 'best' bin, according to some criteria
        for contig, bins in ambiguous_contigs.items():
            # Sort bins by completeness - contamination
            bins.sort(key=lambda bin_id: checkm_stats[bin_id]['completeness'] - checkm_stats[bin_id]['contamination'], reverse=True)
            for bin_id in bins[1:]:  # Skip the first bin since it's considered the best one
                del bin_sets[bin_id][contig]  # Remove contig from bin
    elif keep_ambiguous:
        # If keeping ambiguous, no action needed since they are already listed in multiple bins
        pass

    return bin_sets


# We now need to choose the best version of each bin from all bin refinement iterations
def consolidate_bin_sets(binner_outputs, checkm_stats, output_dir, min_completion, max_contamination, keep_ambiguous, remove_ambiguous):
    # In a dict, for each bin set, create a dict with bins as keys and as value a dict with contig names as keys and contig lengths as values
    bin_sets = get_bin_cont_lengths(binner_outputs)
    
    # Calculate similarity between each pair of bins from each set
    similarity_scores = calculate_weighted_bin_similarity(bin_sets)

    # Convert similarity scores to distances
    distances = {pair: 1 - score for pair, score in similarity_scores.items()}
    # Convert to a square distance matrix
    distance_matrix = squareform(list(distances.values()))

    raise Exception

    # Apply HDBSCAN
    # clusterer = hdbscan.HDBSCAN(min_cluster_size=1, metric='precomputed')
    # clusters = clusterer.fit_predict(distance_matrix)

    # Apply UniDip
    clusterer = UniDip(random_state=0)

    clusters = clusterer.fit(distance_matrix)


    # Relabel clusters labeled -1 (noise) to a unique label
    clusters = [cluster if cluster != -1 else max(clusters) + 1 for cluster in clusters]

    # Build bin to cluster mapping
    bin2cluster = {bin_similarity_tuple[0]: cluster for bin_similarity_tuple, cluster in zip(similarity_scores.keys(), clusters)}

    # Plot clusters
    plot_clusters(distance_matrix, clusters, filename=f"{output_dir}/refined_bin_versions_clusters.svg")

    # Find best bin per cluster
    best_bin_per_cluster = choose_best_bin_per_cluster(bin2cluster, checkm_stats, min_completion, max_contamination)

    # Check if contigs are in more than one bin in the final set
    bin_sets_clean = identify_and_resolve_ambiguous_contigs(bin_sets, best_bin_per_cluster, keep_ambiguous, remove_ambiguous)
    
    return best_bin_per_cluster, bin_sets_clean


def write_consolidated_bins_to_files(bin_sets_updated, best_bin_per_cluster, output_dir, set_name2path):
    for cluster_id, best_bin in best_bin_per_cluster.items():
        bin_set, bin_name = best_bin.split(':')
        # Determine the path for the new consolidated bin file
        output_file = os.path.join(output_dir, f"consolidated_{bin_name}.fasta")

        # Write contigs to the new file
        with open(output_file, 'w') as f_out:
            for contig_id, length in bin_sets_updated[bin_set][bin_name].items():
                # Load the original bin file that contains this contig
                full_parent_path = set_name2path[bin_set]
                org_bin_path = os.path.join(full_parent_path, bin_name)
                id2seq = SeqIO.to_dict(SeqIO.parse(org_bin_path, 'fasta'))
                contig_seq = id2seq[contig_id].seq
                SeqIO.write(SeqIO.SeqRecord(contig_seq, id=contig_id, description=""), f_out, "fasta")


def filter_bins_based_on_checkm2(final_bins_dir, checkm2_stats, min_completion, max_contamination):
    checkm2_results, header = parse_checkm2_output(checkm2_stats, True)
    final_bins = {}
    final_stats_file = os.path.join(final_bins_dir, 'bin_stats.tsv')
    with open(final_stats_file, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for bin_id, metrics in checkm2_results.items():
            if metrics['completeness'] < min_completion or metrics['contamination'] > max_contamination:
                # I am assuing checkm strips the extension and its always fasta
                os.remove(os.path.join(final_bins_dir, f"{bin_id}.fasta"))
            else:
                # Write all bin stats in dictto final stats file
                metrics_str = '\t'.join([str(metrics[header_item]) for header_item in header[1:]])
                f.write(f"{bin_id}\t{metrics_str}\n")
    return final_bins


def main():
    args = parse_arguments()

    # Print parsed arguments to verify (this is just a placeholder, implement your logic here)
    print("Output Directory:", args.output_dir)
    print("Threads:", args.threads)
    print("Min Completion:", args.min_completion)
    print("Max Contamination:", args.max_contamination)
    print("Bin Folders:", args.bin_folders)
    print("Skip Refinement:", args.skip_refinement)
    print("Skip CheckM:", args.skip_checkm)
    print("Skip Consolidation:", args.skip_consolidation)
    print("Keep Ambiguous:", args.keep_ambiguous)
    print("Remove Ambiguous:", args.remove_ambiguous)
    print("Quick:", args.quick)

    # Create symlinks for all combinations of bin folders
    combinations = create_symlinks_and_prepare_combinations(args.bin_folders, args.output_dir)
    
    # Run binning_refiner on each combination in parallel
    with Pool(args.threads) as pool:
        refine_tasks = [(comb, args.output_dir, args.binning_refiner_path) for comb in combinations]
        refined_outputs = pool.starmap(run_binning_refiner, refine_tasks)
    
    # Add original binner outputs to the list of refined outputs
    all_binner_sets = args.bin_folders + refined_outputs

    # After all binning_refiner tasks have completed, run checkm2 iteratively on each refined output
    binning_checkm_stats = run_checkm2_on_refined_outputs(all_binner_sets, args.checkm_db, args.threads)

    best_bin_per_cluster, bin_sets = consolidate_bin_sets(all_binner_sets, binning_checkm_stats, args.output_dir, 
                                                          args.min_completion, args.max_contamination, args.keep_ambiguous, 
                                                          args.remove_ambiguous)

    # Write consolidated bins to file
    final_bins_dir = os.path.join(args.output_dir, 'bins')
    set_name2path = {os.path.basename(binner_folder): binner_folder for binner_folder in all_binner_sets}
    write_consolidated_bins_to_files(bin_sets, best_bin_per_cluster, final_bins_dir, set_name2path)

    # Run CheckM2 on the final set
    run_checkm2(final_bins_dir, args.output_dir, args.checkm_db, args.threads)
        
    # Read CheckM2 output and remove bins that don't meet the criteria
    checkm2_output_file = f"{final_bins_dir}/quality_report.tsv"
    final_bins = filter_bins_based_on_checkm2(final_bins_dir, checkm2_output_file, args.min_completion, args.max_contamination)


if __name__ == "__main__":
    main()
