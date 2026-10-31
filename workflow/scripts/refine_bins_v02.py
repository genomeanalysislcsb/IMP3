#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import itertools
import logging
import shutil
import tempfile
import concurrent.futures
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict, Counter
import hdbscan
from scipy.spatial.distance import squareform, pdist
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import pickle
from pathlib import Path
import glob
import gzip
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class SplitArgs(argparse.Action):
    """Custom action to split comma-separated arguments into lists and expand wildcards."""
    def __call__(self, parser, namespace, values, option_string=None):
        # Split by comma
        paths = values.split(',')
        
        # Expand wildcards in each path
        expanded_paths = []
        for path in paths:
            # Remove any surrounding quotes
            path = path.strip('"\'')
            
            # Check if path contains wildcards
            if '*' in path or '?' in path or '[' in path:
                # Expand wildcards
                glob_results = glob.glob(os.path.expanduser(path))
                if not glob_results:
                    parser.error(f"Wildcard pattern '{path}' did not match any files or directories")
                expanded_paths.extend(glob_results)
            else:
                expanded_paths.append(os.path.expanduser(path))
        
        # Remove duplicates while preserving order
        unique_paths = []
        for path in expanded_paths:
            if path not in unique_paths:
                unique_paths.append(path)
        
        if not unique_paths:
            parser.error(f"No valid paths found in: {values}")
        
        setattr(namespace, self.dest, unique_paths)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Enhanced bin refiner: consolidate and refine metagenomic bins from multiple binning algorithms.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help="Output directory")
    parser.add_argument('-i', '--bin_folders', action=SplitArgs, required=True,
                        help="Comma-separated list of bin folders to use for refinement (supports wildcards, e.g., '/path/to/bins/*')")
    
    # New argument for specifying binner names
    parser.add_argument('--binner-names', type=str, 
                        help="Comma-separated list of names for input bin folders (must match number of folders)")
    
    # Processing parameters
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help="Number of threads")
    parser.add_argument('-b', '--batch_size', type=int, default=None,
                        help="Number of combinations to process in each batch (default: equal to threads for maximum core usage)")
    parser.add_argument('--max_combinations', type=int, default=None,
                        help="Maximum number of bin combinations to process (default: all)")
    
    # Quality thresholds
    parser.add_argument('-c', '--min_completion', type=float, default=70,
                        help="Minimum completion percentage of bins")
    parser.add_argument('-x', '--max_contamination', type=float, default=10,
                        help="Maximum contamination percentage of bins")
    parser.add_argument('--min_bin_size', type=int, default=500,
                        help="Minimum bin size in Kbp for refined bins")
    
    # Contig ambiguity handling
    parser.add_argument('--keep-ambiguous', action='store_true',
                        help="Keep contigs in all bins they appear in")
    parser.add_argument('--remove-ambiguous', action='store_true',
                        help="Remove contigs from all bins if they appear in multiple bins")
    
    # External dependencies
    parser.add_argument('--binning_refiner_path', type=str, required=True,
                        help="Path to binning_refiner script (file path)")
    parser.add_argument('--checkm_db', type=str, required=True,
                        help="Path to CheckM2 database")
    
    # Workflow control
    parser.add_argument('--skip-refinement', action='store_true',
                        help="Skip refinement step (only consolidate original bins)")
    parser.add_argument('--skip-checkm', action='store_true',
                        help="Skip CheckM2 quality assessment")
    parser.add_argument('--skip-consolidation', action='store_true',
                        help="Skip bin consolidation step")
    parser.add_argument('--exclude-singles', action='store_true',
                        help="Exclude single bin sets from refinement (by default, singles are included)")
    
    # Clustering parameters for bin consolidation
    parser.add_argument('--similarity-threshold', type=float, default=0.7,
                        help="Similarity threshold for bin clustering")
    parser.add_argument('--min-cluster-size', type=int, default=2,
                        help="Minimum cluster size for HDBSCAN")
    parser.add_argument('--min-samples', type=int, default=1,
                        help="Minimum samples parameter for HDBSCAN clustering")
    
    # Visualization
    parser.add_argument('--skip-visualization', action='store_true',
                        help="Skip visualization of bin clustering")
    
    # Debug options
    parser.add_argument('--debug', action='store_true',
                        help="Enable debug logging")
    parser.add_argument('--keep-temp', action='store_true',
                        help="Keep temporary files")
    
    args = parser.parse_args()
    
    # Set batch size to threads if not specified
    if args.batch_size is None:
        args.batch_size = args.threads
    
    # Set log level
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    return args


def is_gzipped(file_path):
    """Check if a file is gzipped by looking at its magic bytes."""
    try:
        with open(file_path, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'
    except:
        return False


def open_fasta(filename, mode='r'):
    """Open a FASTA file, automatically handling gzipped files."""
    if is_gzipped(filename):
        return gzip.open(filename, f"{mode}t")  # Open in text mode
    else:
        return open(filename, mode)


def is_fasta_file(filename):
    """Check if a file is a FASTA file based on extension."""
    fasta_extensions = ('.fa', '.fasta', '.fna')
    gz_fasta_extensions = ('.fa.gz', '.fasta.gz', '.fna.gz')
    
    # Check plain FASTA extensions
    if filename.endswith(fasta_extensions):
        return True, False
    
    # Check gzipped FASTA extensions
    if filename.endswith(gz_fasta_extensions):
        return True, True
    
    return False, False


def validate_inputs(args):
    """Validate input parameters and files."""
    # Check if bin folders exist
    logger.info(f"Validating input bin folders: {args.bin_folders}")
    valid_bin_folders = []
    
    for bin_folder in args.bin_folders:
        if not os.path.isdir(bin_folder):
            logger.error(f"Bin folder does not exist: {bin_folder}")
            continue
        
        # Check if bin folder contains bins (including gzipped FASTA)
        bin_files = []
        for filename in os.listdir(bin_folder):
            is_fasta, _ = is_fasta_file(filename)
            if is_fasta:
                bin_files.append(filename)
        
        if not bin_files:
            logger.warning(f"No bin files found in folder: {bin_folder} - skipping")
            continue
        
        valid_bin_folders.append(bin_folder)
    
    if not valid_bin_folders:
        logger.error("No valid bin folders found. Please check paths and try again.")
        sys.exit(1)
    
    # Update args.bin_folders with only valid folders
    args.bin_folders = valid_bin_folders
    logger.info(f"Found {len(valid_bin_folders)} valid bin folders")
    
    # Check binner names if provided
    if args.binner_names:
        binner_names = args.binner_names.split(',')
        if len(binner_names) != len(args.bin_folders):
            logger.error(f"Number of binner names ({len(binner_names)}) must match number of bin folders ({len(args.bin_folders)})")
            sys.exit(1)
    
    # Check for binning_refiner
    if not os.path.isfile(args.binning_refiner_path):
        logger.error(f"binning_refiner not found at: {args.binning_refiner_path}")
        sys.exit(1)
    
    # Check for CheckM2 if not skipped
    if not args.skip_checkm:
        try:
            subprocess.run(["checkm2", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        except FileNotFoundError:
            logger.error("CheckM2 not found in PATH. Please install CheckM2 or use --skip-checkm.")
            sys.exit(1)
        
        # Check if CheckM2 database exists
        if not os.path.exists(args.checkm_db):
            logger.error(f"CheckM2 database not found at: {args.checkm_db}")
            sys.exit(1)
    
    # Check for conflicting ambiguity handling options
    if args.keep_ambiguous and args.remove_ambiguous:
        logger.error("Cannot specify both --keep-ambiguous and --remove-ambiguous")
        sys.exit(1)


def setup_output_directory(output_dir):
    """Create output directory structure."""
    # Create main output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create subdirectories
    subdirs = {
        "combinations": os.path.join(output_dir, "combinations"),
        "refined_bins": os.path.join(output_dir, "refined_bins"),
        "checkm_results": os.path.join(output_dir, "checkm_results"),
        "figures": os.path.join(output_dir, "figures"),
        "final_bins": os.path.join(output_dir, "final_bins"),
        "temp": os.path.join(output_dir, "temp")
    }
    
    for path in subdirs.values():
        os.makedirs(path, exist_ok=True)
    
    return subdirs


def clean_bin_folder_name(folder_path, bin_folder_to_name=None):
    """Clean bin folder names to create better combination names."""
    if bin_folder_to_name and folder_path in bin_folder_to_name:
        return bin_folder_to_name[folder_path]
    
    # Extract basename
    name = os.path.basename(folder_path)
    
    # Remove common suffixes like '_bins', '.bins', '/bins'
    name = re.sub(r'[_\./]bin[s]?$', '', name)
    
    # Remove any other non-alphanumeric characters except underscores
    name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
    
    return name


def generate_bin_combinations(bin_folders, bin_folder_to_name=None, max_combinations=None, min_size=2, exclude_singles=False):
    """Generate combinations of bin folders up to a maximum number."""
    n = len(bin_folders)
    all_combinations = []
    
    # Start with largest combinations and work down
    sizes = list(range(n, min_size - 1, -1))
    if not exclude_singles and min_size > 1:
        sizes.append(1)
    
    # Create sanitized names for bin folders (removing common suffixes)
    sanitized_names = [clean_bin_folder_name(folder, bin_folder_to_name) for folder in bin_folders]
    
    for size in sizes:
        size_combos = []
        for combo_indices in itertools.combinations(range(len(bin_folders)), size):
            combo_paths = [bin_folders[i] for i in combo_indices]
            combo_names = [sanitized_names[i] for i in combo_indices]
            
            size_combos.append({
                'size': size,
                'paths': combo_paths,
                'name': f"combo_size{size}_" + "_".join(combo_names)
            })
        
        # Sort combinations within each size by names for reproducibility
        size_combos.sort(key=lambda x: x['name'])
        all_combinations.extend(size_combos)
    
    # Limit combinations if specified
    if max_combinations is not None and len(all_combinations) > max_combinations:
        logger.warning(f"Limiting to {max_combinations} combinations out of {len(all_combinations)} possible")
        all_combinations = all_combinations[:max_combinations]
    
    return all_combinations


def create_combination_symlinks(combination, output_dir):
    """Create symlinks for a bin combination in the output directory."""
    combo_dir = os.path.join(output_dir, combination['name'])
    os.makedirs(combo_dir, exist_ok=True)
    
    for bin_folder in combination['paths']:
        symlink_name = os.path.join(combo_dir, os.path.basename(bin_folder))
        if not os.path.exists(symlink_name):
            try:
                os.symlink(os.path.abspath(bin_folder), symlink_name)
            except OSError:
                # If symlink fails (e.g., on Windows), copy files instead
                if os.path.exists(symlink_name):
                    shutil.rmtree(symlink_name)
                shutil.copytree(bin_folder, symlink_name)
    
    return combo_dir


def process_bin_combination(combination, output_dir, refiner_script, min_bin_size):
    """Process a single bin combination using binning_refiner."""
    try:
        # Create symlinks for this combination
        combo_dir = create_combination_symlinks(combination, output_dir)
        
        # Prepare output paths
        output_prefix = os.path.join(output_dir, combination['name'] + "_refined")
        
        # Run binning_refiner
        logger.info(f"Running binning_refiner on {combination['name']}")
        
        command = [
            sys.executable,  # Use the current Python interpreter
            refiner_script,  # Direct path to binning_refiner.py
            "-i", combo_dir,
            "-o", output_prefix,
            "-p", combination['name'],
            "-m", str(min_bin_size),
            "-q", "True"  # Silence binning_refiner output
        ]
        
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        # Check if Refined folder was created
        refined_dir = os.path.join(output_prefix, 'Refined')
        if not os.path.exists(refined_dir):
            logger.error(f"Binning_refiner failed to create Refined folder for {combination['name']}")
            return None
        
        # Check if any bins were produced
        bin_files = []
        for filename in os.listdir(refined_dir):
            is_fasta, _ = is_fasta_file(filename)
            if is_fasta:
                bin_files.append(filename)
                
        if not bin_files:
            logger.warning(f"No refined bins produced for {combination['name']}")
            return None
        
        # Rename bin files to include combination name
        for bin_file in bin_files:
            bin_base = os.path.splitext(bin_file)[0]
            bin_ext = os.path.splitext(bin_file)[1]
            
            # Handle .gz extension if present
            if bin_ext == '.gz' and '.' in bin_base:
                bin_base, inner_ext = os.path.splitext(bin_base)
                bin_ext = inner_ext + bin_ext
            
            new_name = f"{combination['name']}_{bin_base}.fasta"
            os.rename(
                os.path.join(refined_dir, bin_file),
                os.path.join(refined_dir, new_name)
            )
        
        logger.info(f"Successfully refined {len(bin_files)} bins for {combination['name']}")
        return refined_dir
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running binning_refiner on {combination['name']}: {e.stderr}")
        return None
    
    except Exception as e:
        logger.error(f"Error processing combination {combination['name']}: {str(e)}")
        return None


def process_combinations_in_parallel(combinations, output_dir, refiner_script, min_bin_size, threads, batch_size):
    """Process bin combinations in parallel batches."""
    all_refined_dirs = []
    total_batches = (len(combinations) + batch_size - 1) // batch_size
    
    # Process combinations in batches
    for i in range(0, len(combinations), batch_size):
        batch = combinations[i:i+batch_size]
        batch_num = i // batch_size + 1
        logger.info(f"Processing batch {batch_num}/{total_batches} with {len(batch)} combinations")
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            future_to_combo = {
                executor.submit(
                    process_bin_combination, 
                    combo, 
                    output_dir, 
                    refiner_script, 
                    min_bin_size
                ): combo for combo in batch
            }
            
            for future in concurrent.futures.as_completed(future_to_combo):
                combo = future_to_combo[future]
                try:
                    refined_dir = future.result()
                    if refined_dir:
                        all_refined_dirs.append(refined_dir)
                except Exception as e:
                    logger.error(f"Exception processing {combo['name']}: {str(e)}")
    
    return all_refined_dirs


def collect_all_bins(bin_dirs):
    """Collect all bin files from multiple directories into a list of paths."""
    all_bins = []
    
    for bin_dir in bin_dirs:
        if not os.path.isdir(bin_dir):
            logger.warning(f"Directory not found: {bin_dir}, skipping")
            continue
            
        bin_files = []
        for filename in os.listdir(bin_dir):
            is_fasta, _ = is_fasta_file(filename)
            if is_fasta:
                bin_files.append(os.path.join(bin_dir, filename))
        
        all_bins.extend(bin_files)
    
    return all_bins


def run_unified_checkm2(all_bin_dirs, output_dir, db_path, threads):
    """Run CheckM2 once on all bins with appropriate naming to track origin."""
    # Create a temporary directory for bin symlinks
    tmp_dir = os.path.join(output_dir, "all_bins_tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Create symlinks to all bins with unique names
    bin_mapping = {}  # Maps new names to original paths
    for bin_dir in all_bin_dirs:
        # Clean up directory name for better identification
        bin_set_name = clean_bin_folder_name(bin_dir)
        
        # Get all bin files (both plain and gzipped FASTA)
        bin_files = []
        for filename in os.listdir(bin_dir):
            is_fasta, _ = is_fasta_file(filename)
            if is_fasta:
                bin_files.append(filename)
        
        for bin_file in bin_files:
            # Create a unique name combining bin set and bin id
            unique_name = f"{bin_set_name}__{bin_file}"
            source_path = os.path.join(bin_dir, bin_file)
            target_path = os.path.join(tmp_dir, unique_name)
            
            # Skip if target file already exists (this can happen with overlapping bin sets)
            if os.path.exists(target_path):
                logger.warning(f"Skipping duplicate bin: {unique_name}")
                continue
                
            # Create symlink (or copy if symlink fails)
            try:
                os.symlink(os.path.abspath(source_path), target_path)
            except OSError as e:
                # Only try to copy if it's not the same file
                if os.path.abspath(source_path) != os.path.abspath(target_path):
                    try:
                        shutil.copy2(source_path, target_path)
                    except shutil.SameFileError:
                        logger.warning(f"Cannot copy {source_path} to {target_path}: same file")
                        continue
                else:
                    logger.warning(f"Cannot symlink {source_path} to {target_path}: {str(e)}")
                    continue
                
            bin_mapping[unique_name] = {
                "bin_set": bin_set_name, 
                "bin_id": bin_file,
                "original_path": os.path.abspath(source_path)
            }
    
    # Run CheckM2 on the temporary directory
    output_folder = os.path.join(output_dir, "unified_checkm2")
    os.makedirs(output_folder, exist_ok=True)
    
    try:
        logger.info(f"Running unified CheckM2 on {len(bin_mapping)} bins")
        
        # Check if db_path points to the database directory or the file
        db_arg = db_path
        if os.path.isdir(db_path):
            # Find the main database file in the directory
            for filename in os.listdir(db_path):
                if filename.endswith('.dmnd') or filename == 'db.dmnd':
                    db_arg = os.path.join(db_path, filename)
                    break
        
        # Build the base command
        command = [
            "checkm2", 
            "predict",
            "--threads", str(threads),
            "--input", tmp_dir,
            "--output-directory", output_folder,
            "--force",
            "--database_path", db_arg
        ]
        
        # Determine file extensions present
        has_fasta = any(f.endswith(('.fasta', '.fa', '.fna')) and not f.endswith('.gz') 
                         for f in os.listdir(tmp_dir))
        has_gz = any(f.endswith('.gz') for f in os.listdir(tmp_dir))
        
        # Add appropriate extension flag based on file types present
        if has_fasta and has_gz:
            logger.info("Detected mixed file extensions (.fasta and .gz)")
            command.extend(["--extension", "fasta,gz"])
        elif has_gz:
            logger.info("Detected gzipped FASTA files only")
            command.extend(["--extension", "gz"])
        elif has_fasta:
            logger.info("Detected FASTA files only")
            command.extend(["--extension", "fasta"])
        else:
            logger.warning("No recognized FASTA files found in temp directory!")
        
        # Log the command for debugging
        logger.debug(f"Running CheckM2 command: {' '.join(command)}")
        
        subprocess.run(command, check=True)
        
        # Save bin mapping for later use
        with open(os.path.join(output_folder, "bin_mapping.pkl"), "wb") as f:
            pickle.dump(bin_mapping, f)
        
        # Parse results and organize by bin set
        checkm_results = {}
        quality_report = os.path.join(output_folder, "quality_report.tsv")
        
        if os.path.exists(quality_report):
            df = pd.read_csv(quality_report, sep='\t')
            
            for _, row in df.iterrows():
                bin_id = row['Bin Id']
                if bin_id in bin_mapping:
                    bin_set = bin_mapping[bin_id]["bin_set"]
                    original_bin_id = bin_mapping[bin_id]["bin_id"]
                    
                    # Initialize bin set in results if not present
                    if bin_set not in checkm_results:
                        checkm_results[bin_set] = {}
                    
                    # Store metrics under original bin ID
                    metrics = {col: row[col] for col in df.columns if col != 'Bin Id'}
                    checkm_results[bin_set][original_bin_id] = metrics
            
            logger.info(f"Parsed CheckM2 results for {sum(len(bins) for bins in checkm_results.values())} bins")
        else:
            logger.error(f"CheckM2 quality report not found at {quality_report}")
        
        return checkm_results, bin_mapping
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running CheckM2: {str(e)}")
        return {}, {}
    
    except Exception as e:
        logger.error(f"Unexpected error running CheckM2: {str(e)}")
        return {}, {}


def extract_contig_info_from_bins(bin_set_path):
    """Extract contig information from all bins in a bin set."""
    bin_contigs = {}
    
    # Get all bin files (both plain and gzipped FASTA)
    bin_files = []
    for filename in os.listdir(bin_set_path):
        is_fasta, _ = is_fasta_file(filename)
        if is_fasta:
            bin_files.append(filename)
    
    for bin_file in bin_files:
        bin_id = bin_file
        bin_contigs[bin_id] = {}
        
        try:
            bin_path = os.path.join(bin_set_path, bin_file)
            with open_fasta(bin_path) as f:
                for contig in SeqIO.parse(f, 'fasta'):
                    contig_id = contig.id
                    contig_length = len(contig.seq)
                    bin_contigs[bin_id][contig_id] = contig_length
        
        except Exception as e:
            logger.error(f"Error processing bin {bin_file}: {str(e)}")
    
    return bin_contigs


def collect_all_bin_contig_info(bin_set_paths, bin_folder_to_name=None):
    """Collect contig information from all bin sets."""
    all_bin_contigs = {}
    
    for bin_set_path in bin_set_paths:
        # Clean up directory name for better identification
        bin_set_name = clean_bin_folder_name(bin_set_path, bin_folder_to_name)
        
        bin_contigs = extract_contig_info_from_bins(bin_set_path)
        
        if bin_contigs:
            all_bin_contigs[bin_set_name] = bin_contigs
        else:
            logger.warning(f"No bin contig information found for {bin_set_name}")
    
    return all_bin_contigs


def calculate_bin_similarities(all_bin_contigs):
    """Calculate similarity between all pairs of bins across all bin sets."""
    # Create a flat list of all bins with their contigs
    all_bins = []
    bin_labels = []
    
    for bin_set, bin_contigs in all_bin_contigs.items():
        for bin_id, contigs in bin_contigs.items():
            all_bins.append((bin_set, bin_id, contigs))
            bin_labels.append(f"{bin_set}:{bin_id}")
    
    n_bins = len(all_bins)
    similarity_matrix = np.zeros((n_bins, n_bins))
    
    # Calculate pairwise similarity using a more comprehensive metric
    for i in range(n_bins):
        bin_set_i, bin_id_i, contigs_i = all_bins[i]
        bin_i_contigs = set(contigs_i.keys())
        bin_i_length = sum(contigs_i.values())
        
        # Set diagonal to 1 (each bin is identical to itself)
        similarity_matrix[i, i] = 1.0
        
        for j in range(i+1, n_bins):
            bin_set_j, bin_id_j, contigs_j = all_bins[j]
            bin_j_contigs = set(contigs_j.keys())
            bin_j_length = sum(contigs_j.values())
            
            # Skip comparison if bins are from the same set (never cluster these)
            if bin_set_i == bin_set_j:
                similarity_matrix[i, j] = 0.0
                similarity_matrix[j, i] = 0.0
                continue
            
            # Calculate shared contigs and their lengths
            shared_contigs = bin_i_contigs.intersection(bin_j_contigs)
            
            if not shared_contigs:
                similarity_matrix[i, j] = 0.0
                similarity_matrix[j, i] = 0.0
                continue
            
            # Calculate total shared length and Jaccard similarity weighted by length
            shared_length = sum(min(contigs_i[contig], contigs_j[contig]) for contig in shared_contigs)
            
            # Similarity is 2 * shared length / total length
            if bin_i_length + bin_j_length > 0:
                similarity = 2 * shared_length / (bin_i_length + bin_j_length)
            else:
                similarity = 0
            
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
    
    return similarity_matrix, bin_labels


def cluster_bins(similarity_matrix, bin_labels, min_cluster_size=2, min_samples=1, cluster_selection_epsilon=0.05):
    """Cluster bins based on similarity using HDBSCAN."""
    # Ensure min_cluster_size is at least 2
    min_cluster_size = max(2, min_cluster_size)
    
    # Handle the case where there are too few samples
    if len(bin_labels) < min_cluster_size:
        logger.warning(f"Only {len(bin_labels)} bins available, not enough for clustering with min_cluster_size={min_cluster_size}")
        labels = np.array([-i-1 for i in range(len(bin_labels))])  # Assign each bin to its own cluster
    else:
        # Convert similarities to distances
        distance_matrix = 1 - similarity_matrix
        
        # Apply HDBSCAN clustering
        try:
            clusterer = hdbscan.HDBSCAN(
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                metric='precomputed',
                cluster_selection_epsilon=cluster_selection_epsilon,
                cluster_selection_method='eom'  # Excess of Mass for more conservative clustering
            )
            
            labels = clusterer.fit_predict(distance_matrix)
        except ValueError as e:
            if "Min cluster size must be greater than one" in str(e):
                logger.warning("HDBSCAN error: Min cluster size must be greater than one. Creating singleton clusters.")
                labels = np.array([-i-1 for i in range(len(bin_labels))])  # Assign each bin to its own cluster
            else:
                raise  # Re-raise if it's a different ValueError
    
    # Map bins to clusters
    bin_to_cluster = {bin_label: cluster_id for bin_label, cluster_id in zip(bin_labels, labels)}
    
    # Group bins by cluster
    clusters = defaultdict(list)
    for bin_label, cluster_id in bin_to_cluster.items():
        if cluster_id != -1:  # Skip noise points for now
            clusters[cluster_id].append(bin_label)
    
    # Handle unclustered bins (noise points)
    noise_points = [bin_label for bin_label, cluster_id in bin_to_cluster.items() if cluster_id == -1]
    
    # Assign each noise point to a singleton cluster
    for i, bin_label in enumerate(noise_points):
        unique_id = -i - 1000  # Create unique negative IDs starting from -1000
        bin_to_cluster[bin_label] = unique_id
        clusters[unique_id] = [bin_label]
    
    logger.info(f"Found {len(clusters) - len(noise_points)} non-singleton clusters and {len(noise_points)} singleton clusters")
    
    return bin_to_cluster, clusters


def validate_clusters(clusters, all_bin_contigs, similarity_matrix, bin_labels, min_similarity=0.4):
    """Validate clusters to ensure they are coherent."""
    validated_clusters = {}
    bin_indices = {label: idx for idx, label in enumerate(bin_labels)}
    
    for cluster_id, bin_labels_in_cluster in clusters.items():
        # Skip singleton clusters (they're already valid)
        if len(bin_labels_in_cluster) <= 1:
            validated_clusters[cluster_id] = bin_labels_in_cluster
            continue
        
        # Check if all bins in cluster are sufficiently similar to each other
        is_valid = True
        for i, bin_i in enumerate(bin_labels_in_cluster):
            idx_i = bin_indices[bin_i]
            
            for bin_j in bin_labels_in_cluster[i+1:]:
                idx_j = bin_indices[bin_j]
                similarity = similarity_matrix[idx_i, idx_j]
                
                if similarity < min_similarity:
                    is_valid = False
                    logger.warning(f"Cluster {cluster_id} contains dissimilar bins: {bin_i} and {bin_j} with similarity {similarity:.3f}")
                    break
            
            if not is_valid:
                break
        
        if is_valid:
            validated_clusters[cluster_id] = bin_labels_in_cluster
        else:
            # Split invalid cluster into singleton clusters
            for i, bin_label in enumerate(bin_labels_in_cluster):
                unique_id = -(1000000 + cluster_id*1000 + i)  # Create deterministic unique negative IDs
                validated_clusters[unique_id] = [bin_label]
                logger.info(f"Split bin {bin_label} from invalid cluster {cluster_id} into singleton cluster {unique_id}")
    
    return validated_clusters


def visualize_bin_clustering(similarity_matrix, bin_labels, bin_to_cluster, output_path):
    """Visualize bin clustering using UMAP."""
    try:
        # Convert similarity matrix to distance matrix
        distance_matrix = 1 - similarity_matrix
        
        # Apply UMAP for dimensionality reduction
        reducer = umap.UMAP(
            n_neighbors=min(15, len(bin_labels) - 1),
            min_dist=0.1,
            metric='precomputed',
            random_state=42  # For reproducibility
        )
        embedding = reducer.fit_transform(distance_matrix)
        
        # Get cluster labels for each bin
        clusters = [bin_to_cluster[label] for label in bin_labels]
        
        # Extract bin set and bin id for each point
        bin_sets = []
        bin_ids = []
        for label in bin_labels:
            parts = label.split(':')
            bin_sets.append(parts[0])
            if len(parts) > 1:
                bin_ids.append(parts[1])
            else:
                bin_ids.append("unknown")
        
        # Create a DataFrame for easier plotting
        df = pd.DataFrame({
            'UMAP1': embedding[:, 0],
            'UMAP2': embedding[:, 1],
            'Cluster': clusters,
            'BinSet': bin_sets,
            'BinID': bin_ids,
            'Label': bin_labels
        })
        
        # Create visualization
        plt.figure(figsize=(14, 12))
        
        # Create a scatter plot with bin sets shown by different markers
        unique_bin_sets = sorted(set(bin_sets))
        markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', '|', '_']
        
        # Plot points by cluster
        scatter = plt.scatter(
            df['UMAP1'],
            df['UMAP2'],
            c=df['Cluster'],
            cmap='tab20',
            s=100,
            alpha=0.7,
            edgecolors='k',
            linewidths=0.5
        )
        
        # Add labels for clusters with more than one bin
        cluster_counts = Counter(clusters)
        for cluster, count in cluster_counts.items():
            if count > 1:
                points = df[df['Cluster'] == cluster]
                center_x = points['UMAP1'].mean()
                center_y = points['UMAP2'].mean()
                plt.annotate(
                    f"C{cluster}",
                    (center_x, center_y),
                    fontsize=12,
                    fontweight='bold',
                    ha='center',
                    va='center',
                    bbox=dict(boxstyle="circle,pad=0.3", fc="white", ec="black", alpha=0.8)
                )
        
        # Add a legend for bin sets
        for i, bin_set in enumerate(unique_bin_sets):
            marker = markers[i % len(markers)]
            plt.scatter(
                [], [], 
                marker=marker, 
                s=100, 
                c='gray', 
                alpha=0.7,
                edgecolors='k',
                linewidths=0.5,
                label=bin_set
            )
        
        plt.legend(title="Bin Sets", loc="upper right")
        
        # Add interactive hover labels (bin_id) for points
        for i, row in df.iterrows():
            plt.annotate(
                row['BinID'],
                (row['UMAP1'], row['UMAP2']),
                fontsize=8,
                alpha=0,  # Make text invisible initially
                textcoords="offset points",
                xytext=(0, 5),
                ha='center'
            )
        
        plt.title('Bin Clustering Visualization')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        # Also save an interactive HTML version using plotly if available
        try:
            import plotly.express as px
            import plotly.graph_objects as go
            
            fig = px.scatter(
                df, x='UMAP1', y='UMAP2', 
                color='Cluster', symbol='BinSet',
                hover_data=['BinID', 'BinSet', 'Cluster'],
                title='Bin Clustering Visualization (Interactive)'
            )
            
            # Add cluster labels
            for cluster, count in cluster_counts.items():
                if count > 1:
                    points = df[df['Cluster'] == cluster]
                    center_x = points['UMAP1'].mean()
                    center_y = points['UMAP2'].mean()
                    fig.add_annotation(
                        x=center_x, y=center_y,
                        text=f"C{cluster}",
                        showarrow=False,
                        font=dict(size=14, color="black"),
                        bgcolor="white",
                        bordercolor="black",
                        borderwidth=1
                    )
            
            html_path = output_path.replace('.png', '.html')
            fig.write_html(html_path)
            logger.info(f"Saved interactive bin clustering visualization to {html_path}")
            
        except ImportError:
            logger.debug("plotly not available, skipping interactive visualization")
        
        plt.close()
        logger.info(f"Saved bin clustering visualization to {output_path}")
        return True
    
    except Exception as e:
        logger.error(f"Error creating bin clustering visualization: {str(e)}")
        return False


def select_best_bins(bin_to_cluster, clusters, checkm_results, min_completion, max_contamination):
    """Select the best bin from each cluster based on CheckM2 results."""
    best_bins = {}
    
    for cluster_id, bin_labels in clusters.items():
        best_bin = None
        best_score = -float('inf')
        
        for bin_label in bin_labels:
            bin_set, bin_id = bin_label.split(':')
            
            # Skip if CheckM2 results not available for this bin
            if bin_set not in checkm_results or bin_id not in checkm_results[bin_set]:
                logger.warning(f"CheckM2 results not available for {bin_label}")
                continue
            
            # Get bin quality metrics
            metrics = checkm_results[bin_set][bin_id]
            try:
                completeness = float(metrics['completeness'])
                contamination = float(metrics['contamination'])
            except (KeyError, ValueError) as e:
                logger.warning(f"Cannot get completeness/contamination for {bin_label}: {str(e)}")
                continue
            
            # Check if bin meets quality thresholds
            if completeness < min_completion or contamination > max_contamination:
                logger.debug(f"Bin {bin_label} doesn't meet quality thresholds: {completeness:.1f}% completeness, {contamination:.1f}% contamination")
                continue
            
            # Calculate score: completeness - 5*contamination
            # This weighting penalizes contamination more heavily
            score = completeness - 5 * contamination
            
            if score > best_score:
                best_score = score
                best_bin = bin_label
        
        if best_bin:
            best_bins[cluster_id] = best_bin
            logger.debug(f"Selected {best_bin} as best bin for cluster {cluster_id} with score {best_score:.1f}")
        else:
            logger.warning(f"No bin in cluster {cluster_id} meets quality thresholds")
    
    logger.info(f"Selected {len(best_bins)} best bins from {len(clusters)} clusters")
    return best_bins


def resolve_ambiguous_contigs(best_bins, all_bin_contigs, checkm_results, keep_ambiguous=False, remove_ambiguous=False):
    """Resolve contigs that appear in multiple bins.
    
    By default, ambiguous contigs are kept only in the bin with the highest quality score.
    If keep_ambiguous=True, contigs remain in all bins.
    If remove_ambiguous=True, ambiguous contigs are removed from all bins.
    """
    # Extract contig information for best bins
    best_bin_contigs = {}
    for cluster_id, bin_label in best_bins.items():
        bin_set, bin_id = bin_label.split(':')
        if bin_set in all_bin_contigs and bin_id in all_bin_contigs[bin_set]:
            best_bin_contigs[bin_label] = all_bin_contigs[bin_set][bin_id]
        else:
            logger.warning(f"Contig information not found for {bin_label}")
    
    # Identify ambiguous contigs
    contig_to_bins = defaultdict(list)
    for bin_label, contigs in best_bin_contigs.items():
        for contig_id in contigs:
            contig_to_bins[contig_id].append(bin_label)
    
    ambiguous_contigs = {contig: bins for contig, bins in contig_to_bins.items() if len(bins) > 1}
    logger.info(f"Found {len(ambiguous_contigs)} ambiguous contigs in {len(best_bin_contigs)} bins")
    
    # Create a new dictionary with resolved contigs
    resolved_bin_contigs = {bin_label: contigs.copy() for bin_label, contigs in best_bin_contigs.items()}
    
    # Resolve ambiguity based on strategy
    if remove_ambiguous:
        # Remove ambiguous contigs from all bins
        for contig_id, bin_labels in ambiguous_contigs.items():
            for bin_label in bin_labels:
                if bin_label in resolved_bin_contigs and contig_id in resolved_bin_contigs[bin_label]:
                    del resolved_bin_contigs[bin_label][contig_id]
        
        logger.info(f"Removed all {len(ambiguous_contigs)} ambiguous contigs from all bins")
    
    elif not keep_ambiguous:
        # Default: keep ambiguous contigs only in the best bin (highest quality score)
        for contig_id, bin_labels in ambiguous_contigs.items():
            # Calculate quality score for each bin
            bin_scores = {}
            for bin_label in bin_labels:
                bin_set, bin_id = bin_label.split(':')
                if bin_set in checkm_results and bin_id in checkm_results[bin_set]:
                    metrics = checkm_results[bin_set][bin_id]
                    try:
                        completeness = float(metrics['completeness'])
                        contamination = float(metrics['contamination'])
                        # Quality score: completeness - 5*contamination
                        score = completeness - 5 * contamination
                        bin_scores[bin_label] = score
                    except (KeyError, ValueError):
                        bin_scores[bin_label] = -float('inf')
                else:
                    bin_scores[bin_label] = -float('inf')
            
            # Find the bin with the highest score
            if bin_scores:
                best_bin = max(bin_scores.items(), key=lambda x: x[1])[0]
                
                # Remove the contig from all other bins
                for bin_label in bin_labels:
                    if bin_label != best_bin and bin_label in resolved_bin_contigs and contig_id in resolved_bin_contigs[bin_label]:
                        del resolved_bin_contigs[bin_label][contig_id]
        
        logger.info(f"Kept {len(ambiguous_contigs)} ambiguous contigs only in their best bins")
    
    # If keep_ambiguous, no action needed - contigs remain in all bins
    if keep_ambiguous:
        logger.info(f"Kept all {len(ambiguous_contigs)} ambiguous contigs in all bins")
    
    return resolved_bin_contigs


def create_final_bins(resolved_bin_contigs, bin_mapping, output_dir):
    """Create final bin files from resolved contig information."""
    os.makedirs(output_dir, exist_ok=True)
    
    bin_count = 0
    bin_info = {}
    
    for bin_label, contigs in resolved_bin_contigs.items():
        if not contigs:
            logger.warning(f"Bin {bin_label} has no contigs after ambiguity resolution, skipping")
            continue
        
        bin_set, bin_id = bin_label.split(':')
        output_bin_path = os.path.join(output_dir, f"refined_bin_{bin_count:03d}.fasta")
        
        try:
            # Find the original bin path using the bin mapping
            original_path = None
            for map_key, map_info in bin_mapping.items():
                if map_info["bin_set"] == bin_set and map_info["bin_id"] == bin_id:
                    original_path = map_info["original_path"]
                    break
            
            if not original_path:
                logger.error(f"Could not find original path for bin {bin_label}")
                continue
            
            # Create the new bin file with only the selected contigs
            with open(output_bin_path, 'w') as out_file:
                with open_fasta(original_path) as in_file:
                    for record in SeqIO.parse(in_file, 'fasta'):
                        if record.id in contigs:
                            SeqIO.write(record, out_file, 'fasta')
            
            # Store bin information
            bin_info[f"refined_bin_{bin_count:03d}"] = {
                "original_bin": bin_label,
                "contig_count": len(contigs),
                "size_bp": sum(contigs.values())
            }
            
            bin_count += 1
            logger.debug(f"Created bin file: {output_bin_path} with {len(contigs)} contigs")
        
        except Exception as e:
            logger.error(f"Error creating final bin for {bin_label}: {str(e)}")
    
    logger.info(f"Created {bin_count} final bins in {output_dir}")
    
    # Save bin information to file
    with open(os.path.join(output_dir, "bin_info.pkl"), "wb") as f:
        pickle.dump(bin_info, f)
    
    return bin_count, bin_info


def create_final_report(best_bins, checkm_results, bin_info, output_path):
    """Create a final report with bin quality metrics."""
    try:
        with open(output_path, 'w') as f:
            f.write("Bin_ID\tOriginal_Bin\tCompleteness\tContamination\tGenome_Size_bp\tContig_Count\n")
            
            for bin_id, info in bin_info.items():
                original_bin = info["original_bin"]
                bin_set, bin_name = original_bin.split(':')
                
                if bin_set in checkm_results and bin_name in checkm_results[bin_set]:
                    metrics = checkm_results[bin_set][bin_name]
                    completeness = metrics.get('completeness', 'N/A')
                    contamination = metrics.get('contamination', 'N/A')
                    
                    # Use the stored size in bin_info for consistency
                    genome_size = info["size_bp"]
                    contig_count = info["contig_count"]
                    
                    f.write(f"{bin_id}\t{original_bin}\t{completeness}\t{contamination}\t{genome_size}\t{contig_count}\n")
        
        logger.info(f"Created final report at {output_path}")
        return True
    
    except Exception as e:
        logger.error(f"Error creating final report: {str(e)}")
        return False


def plot_bin_quality(best_bins, checkm_results, bin_info, output_path):
    """Create a plot of bin quality metrics."""
    try:
        data = []
        
        for bin_id, info in bin_info.items():
            original_bin = info["original_bin"]
            bin_set, bin_name = original_bin.split(':')
            
            if bin_set in checkm_results and bin_name in checkm_results[bin_set]:
                metrics = checkm_results[bin_set][bin_name]
                try:
                    completeness = float(metrics['completeness'])
                    contamination = float(metrics['contamination'])
                    genome_size_mb = info["size_bp"] / 1_000_000  # Convert to Mb
                    
                    data.append({
                        'Bin': bin_id,
                        'Completeness': completeness,
                        'Contamination': contamination,
                        'Genome_Size_Mb': genome_size_mb,
                        'Original_Bin': original_bin
                    })
                except (KeyError, ValueError):
                    logger.warning(f"Missing or invalid metrics for {bin_id}")
        
        if not data:
            logger.warning("No valid bin quality data to plot")
            return False
        
        # Create a DataFrame
        df = pd.DataFrame(data)
        
        # Sort by completeness
        df = df.sort_values('Completeness', ascending=False)
        
        # Create the plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [3, 1]})
        
        # Set up bar positions
        x = np.arange(len(df))
        width = 0.35
        
        # Plot completeness and contamination
        bars1 = ax1.bar(x - width/2, df['Completeness'], width, label='Completeness (%)', color='#3274A1')
        bars2 = ax1.bar(x + width/2, df['Contamination'], width, label='Contamination (%)', color='#E1812C')
        
        # Add value labels
        for bar in bars1:
            height = bar.get_height()
            ax1.annotate(f'{height:.1f}',
                         xy=(bar.get_x() + bar.get_width() / 2, height),
                         xytext=(0, 3),  # 3 points vertical offset
                         textcoords="offset points",
                         ha='center', va='bottom',
                         fontsize=8)
        
        for bar in bars2:
            height = bar.get_height()
            if height > 1:  # Only show contamination if it's > 1%
                ax1.annotate(f'{height:.1f}',
                             xy=(bar.get_x() + bar.get_width() / 2, height),
                             xytext=(0, 3),  # 3 points vertical offset
                             textcoords="offset points",
                             ha='center', va='bottom',
                             fontsize=8)
        
        # Plot genome size on the second axis
        ax2.bar(x, df['Genome_Size_Mb'], width, label='Genome Size (Mb)', color='#7570B3')
        
        # Set labels and title
        ax1.set_ylabel('Percentage (%)')
        ax1.set_title('Bin Quality Assessment')
        ax1.set_ylim(0, 105)  # Allow for annotation space
        ax1.legend()
        
        ax2.set_xlabel('Bin')
        ax2.set_ylabel('Size (Mb)')
        
        # Set x-ticks at the center of the groups
        ax1.set_xticks(x)
        ax1.set_xticklabels([])  # No labels on the first subplot
        ax2.set_xticks(x)
        ax2.set_xticklabels(df['Bin'], rotation=90)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Created bin quality plot at {output_path}")
        return True
    
    except Exception as e:
        logger.error(f"Error creating bin quality plot: {str(e)}")
        return False


def run_checkm2_on_final_bins(input_dir, output_dir, db_path, threads):
    """Run CheckM2 on the final bin set."""
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        logger.info(f"Running CheckM2 on final bins in {input_dir}")
        
        # Check if db_path points to the database directory or the file
        db_arg = db_path
        if os.path.isdir(db_path):
            # Find the main database file in the directory
            for filename in os.listdir(db_path):
                if filename.endswith('.dmnd') or filename == 'db.dmnd':
                    db_arg = os.path.join(db_path, filename)
                    break
        
        command = [
            "checkm2", 
            "predict",
            "--threads", str(threads),
            "--input", input_dir,
            "--output-directory", output_dir,
            "--extension", "fasta",
            "--force",
            "--database_path", db_arg
        ]
        
        subprocess.run(command, check=True)
        
        # Check if quality report was created
        quality_report = os.path.join(output_dir, "quality_report.tsv")
        if not os.path.exists(quality_report):
            logger.error(f"CheckM2 did not produce a quality report for final bins")
            return False
        
        logger.info(f"Completed CheckM2 analysis for final bins")
        return True
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running CheckM2 on final bins: {str(e)}")
        return False
    
    except Exception as e:
        logger.error(f"Unexpected error running CheckM2 on final bins: {str(e)}")
        return False


def main():
    """Main function coordinating the workflow."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Validate inputs
    validate_inputs(args)
    
    # Create the main output directory first
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create binner names mapping
    if args.binner_names:
        binner_names = args.binner_names.split(',')
        if len(binner_names) != len(args.bin_folders):
            logger.error(f"Number of binner names ({len(binner_names)}) must match number of bin folders ({len(args.bin_folders)})")
            sys.exit(1)
    else:
        # Generate sequential names
        binner_names = [f"binner_{i+1:03d}" for i in range(len(args.bin_folders))]
    
    # Create a mapping between input paths and binner names
    bin_folder_to_name = {folder: name for folder, name in zip(args.bin_folders, binner_names)}
    
    # Save the mapping for reference (now that the directory exists)
    mapping_file = os.path.join(args.output_dir, "binner_mapping.tsv")
    try:
        with open(mapping_file, 'w') as f:
            f.write("Bin_Folder\tBinner_Name\n")
            for folder, name in bin_folder_to_name.items():
                f.write(f"{folder}\t{name}\n")
        logger.info(f"Saved binner mapping to {mapping_file}")
    except Exception as e:
        logger.error(f"Error writing mapping file {mapping_file}: {str(e)}")
    
    # Set up output directory structure
    output_dirs = setup_output_directory(args.output_dir)
    
    # Create a temporary file to indicate keeping temp files if requested
    if args.keep_temp:
        Path(os.path.join(args.output_dir, "keep_temp")).touch()
    
    # Generate bin combinations
    logger.info(f"Generating combinations from {len(args.bin_folders)} bin folders")
    combinations = generate_bin_combinations(
        args.bin_folders,
        bin_folder_to_name=bin_folder_to_name,
        max_combinations=args.max_combinations,
        min_size=2,
        exclude_singles=args.exclude_singles
    )
    logger.info(f"Generated {len(combinations)} combinations")
    
    # Process combinations
    all_bin_dirs = []
    if not args.skip_refinement:
        logger.info("Starting refinement process")
        refined_dirs = process_combinations_in_parallel(
            combinations,
            output_dirs["combinations"],
            args.binning_refiner_path,
            args.min_bin_size,
            args.threads,
            args.batch_size
        )
        all_bin_dirs.extend(refined_dirs)
    
    # Add original bin folders
    all_bin_dirs.extend(args.bin_folders)
    logger.info(f"Total bin directories to process: {len(all_bin_dirs)}")
    
    # Run CheckM2
    checkm_results = {}
    bin_mapping = {}
    if not args.skip_checkm:
        logger.info("Running CheckM2 on all bin sets")
        checkm_results, bin_mapping = run_unified_checkm2(
            all_bin_dirs,
            output_dirs["checkm_results"],
            args.checkm_db,
            args.threads
        )
    
    # Consolidate bins
    if not args.skip_consolidation:
        logger.info("Consolidating bins")
        
        # Collect contig information
        all_bin_contigs = collect_all_bin_contig_info(all_bin_dirs, bin_folder_to_name)
        logger.info(f"Collected contig information for {len(all_bin_contigs)} bin sets")
        
        # Calculate bin similarities
        logger.info("Calculating bin similarities")
        similarity_matrix, bin_labels = calculate_bin_similarities(all_bin_contigs)
        
        # Cluster bins
        logger.info("Clustering bins")
        bin_to_cluster, clusters = cluster_bins(
            similarity_matrix,
            bin_labels,
            min_cluster_size=args.min_cluster_size,
            min_samples=args.min_samples,
            cluster_selection_epsilon=0.05  # Conservative clustering
        )
        
        # Validate clusters
        logger.info("Validating clusters")
        validated_clusters = validate_clusters(
            clusters,
            all_bin_contigs,
            similarity_matrix,
            bin_labels,
            min_similarity=0.4
        )
        
        # Visualize clustering
        if not args.skip_visualization:
            logger.info("Creating bin clustering visualization")
            visualize_bin_clustering(
                similarity_matrix,
                bin_labels,
                bin_to_cluster,
                os.path.join(output_dirs["figures"], "bin_clustering.png")
            )
        
        # Select best bins
        logger.info("Selecting best bins from each cluster")
        best_bins = select_best_bins(
            bin_to_cluster,
            validated_clusters,
            checkm_results,
            args.min_completion,
            args.max_contamination
        )
        
        # Resolve ambiguous contigs
        logger.info("Resolving ambiguous contigs")
        resolved_bin_contigs = resolve_ambiguous_contigs(
            best_bins,
            all_bin_contigs,
            checkm_results,
            keep_ambiguous=args.keep_ambiguous,
            remove_ambiguous=args.remove_ambiguous
        )
        
        # Create final bins
        logger.info("Creating final bin files")
        bin_count, bin_info = create_final_bins(
            resolved_bin_contigs,
            bin_mapping,
            output_dirs["final_bins"]
        )
        
        # Create final report
        logger.info("Creating final report")
        create_final_report(
            best_bins,
            checkm_results,
            bin_info,
            os.path.join(args.output_dir, "final_bins_report.tsv")
        )
        
        # Plot bin quality
        logger.info("Creating bin quality plot")
        plot_bin_quality(
            best_bins,
            checkm_results,
            bin_info,
            os.path.join(output_dirs["figures"], "bin_quality.png")
        )
        
        # Run CheckM2 on final bins
        if not args.skip_checkm and bin_count > 0:
            logger.info("Running CheckM2 on final bins")
            run_checkm2_on_final_bins(
                output_dirs["final_bins"],
                os.path.join(args.output_dir, "final_bins_checkm2"),
                args.checkm_db,
                args.threads
            )
    
    # Clean up temporary files
    if not args.keep_temp:
        logger.info("Cleaning up temporary files")
        try:
            shutil.rmtree(output_dirs["temp"])
            temp_files = [os.path.join(args.output_dir, "keep_temp")]
            for file in temp_files:
                if os.path.exists(file):
                    os.remove(file)
        except Exception as e:
            logger.warning(f"Error cleaning up temporary files: {str(e)}")
    
    logger.info("Bin refinement completed successfully")
    

if __name__ == "__main__":
    main()