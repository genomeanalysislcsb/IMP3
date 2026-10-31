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
import hashlib # Added for filename hashing

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Define filename length limit (common limit, adjust if needed for specific filesystem)
FILENAME_MAX_LEN = 255


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

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MODIFIED FUNCTION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
                        help="Comma-separated list of input bin folders (supports wildcards, e.g., '/path/to/bins/*')")

    # New argument for specifying binner names
    parser.add_argument('--binner-names', type=str,
                        help="Comma-separated list of names for input bin folders (must match number of folders)")

    # Processing parameters
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help="Number of threads")
    parser.add_argument('-b', '--batch_size', type=int, default=None,
                        help="Number of combinations to process in each batch (default: equal to threads for maximum core usage)")
    parser.add_argument('--max_combinations', type=int, default=None,
                        help="Maximum number of bin combinations (size >= 2) to process (default: all)")

    # Quality thresholds (applied during consolidation)
    parser.add_argument('-c', '--min_completion', type=float, default=70,
                        help="Minimum completion percentage for final bins")
    parser.add_argument('-x', '--max_contamination', type=float, default=10,
                        help="Maximum contamination percentage for final bins")
    parser.add_argument('--min_bin_size', type=int, default=500,
                        help="Minimum bin size in Kbp for final bins")

    # Contig ambiguity handling (applied during consolidation)
    parser.add_argument('--keep-ambiguous', action='store_true',
                        help="Keep contigs in all final bins they appear in")
    parser.add_argument('--remove-ambiguous', action='store_true',
                        help="Remove contigs from all final bins if they appear in multiple bins")

    # External dependencies
    parser.add_argument('--binning_refiner_path', type=str, required=True,
                        help="Path to binning_refiner script (e.g., binning_refiner_v002.py)")
    parser.add_argument('--checkm_db', type=str, required=True,
                        help="Path to CheckM2 database")

    # Workflow control
    parser.add_argument('--skip-refinement', action='store_true',
                        help="Skip refinement step (no bin combinations processed)")
    parser.add_argument('--skip-checkm', action='store_true',
                        help="Skip CheckM2 quality assessment (runs on original + refined bins if refinement not skipped)")
    parser.add_argument('--skip-consolidation', action='store_true',
                        help="Skip bin consolidation step (clustering, selection, ambiguity handling)")
    # Removed --exclude-singles as combinations always start at size 2 now

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
                        help="Keep temporary files (including refinement outputs)")

    args = parser.parse_args()

    # Set batch size to threads if not specified
    if args.batch_size is None:
        args.batch_size = args.threads

    # Set log level
    if args.debug:
        logger.setLevel(logging.DEBUG)

    # Ensure min_cluster_size is valid
    if args.min_cluster_size <= 1:
         parser.error("--min-cluster-size must be > 1")


    return args
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END MODIFIED FUNCTION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def is_gzipped(file_path):
    """Check if a file is gzipped by looking at its magic bytes."""
    try:
        with open(file_path, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'
    except:
        return False


def open_fasta(filename, mode='r', encoding='utf-8', errors='ignore'):
    """Open a FASTA file, automatically handling gzipped files and encoding."""
    if is_gzipped(filename):
        # Use gzip.open in text mode with specified encoding and error handling
        return gzip.open(filename, f"{mode}t", encoding=encoding, errors=errors)
    else:
        # Use standard open with specified encoding and error handling
        return open(filename, mode, encoding=encoding, errors=errors)


def is_fasta_file(filename):
    """Check if a file is a FASTA file based on extension."""
    # Convert to lowercase for case-insensitive matching
    lower_filename = filename.lower()
    fasta_extensions = ('.fa', '.fasta', '.fna')
    gz_fasta_extensions = ('.fa.gz', '.fasta.gz', '.fna.gz')

    # Check gzipped FASTA extensions first
    if any(lower_filename.endswith(ext) for ext in gz_fasta_extensions):
        return True, True

    # Check plain FASTA extensions
    if any(lower_filename.endswith(ext) for ext in fasta_extensions):
        return True, False

    return False, False


def validate_inputs(args):
    """Validate input parameters and files."""
    # Check if bin folders exist
    logger.info(f"Validating input bin folders: {args.bin_folders}")
    valid_bin_folders = []

    for bin_folder in args.bin_folders:
        bin_folder = os.path.abspath(os.path.expanduser(bin_folder)) # Ensure absolute path
        if not os.path.isdir(bin_folder):
            logger.error(f"Bin folder does not exist: {bin_folder}")
            continue


        # Check if bin folder contains bins (including gzipped FASTA)
        bin_files = []
        try:
             for filename in os.listdir(bin_folder):
                 is_fasta, _ = is_fasta_file(filename)
                 if is_fasta:
                     bin_files.append(filename)
        except OSError as e:
            logger.error(f"Cannot access bin folder content: {bin_folder}. Error: {e}")
            continue # Skip this folder

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
    if not args.skip_refinement: # Only required if refinement is not skipped
        if not os.path.isfile(args.binning_refiner_path):
            logger.error(f"binning_refiner not found at: {args.binning_refiner_path}")
            sys.exit(1)
        if not os.access(args.binning_refiner_path, os.X_OK):
             logger.warning(f"binning_refiner script at {args.binning_refiner_path} may not be executable.")


    # Check for CheckM2 if not skipped
    if not args.skip_checkm:
        try:
            # Check if checkm2 executable is found
            checkm2_path = shutil.which("checkm2")
            if not checkm2_path:
                raise FileNotFoundError("checkm2 executable not found in PATH")
            # Run checkm2 --version to ensure it's working
            logger.debug(f"Found checkm2 executable at: {checkm2_path}")
            subprocess.run([checkm2_path, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except (FileNotFoundError, subprocess.CalledProcessError) as e:
            logger.error(f"CheckM2 not found in PATH or failed to execute. Please install CheckM2 or use --skip-checkm. Error: {e}")
            sys.exit(1)

        # Check if CheckM2 database exists
        if not os.path.exists(args.checkm_db):
            logger.error(f"CheckM2 database not found at: {args.checkm_db}")
            sys.exit(1)

    # Check for conflicting ambiguity handling options
    if args.keep_ambiguous and args.remove_ambiguous:
        logger.error("Cannot specify both --keep-ambiguous and --remove-ambiguous")
        sys.exit(1)


def setup_output_directory(output_dir, args): # Added 'args' parameter
    """Create output directory structure."""
    # Create main output directory
    os.makedirs(output_dir, exist_ok=True)

    # Create subdirectories
    subdirs = {
        #"combinations": os.path.join(output_dir, "combinations"), # No longer needed? Temp holds inputs.
        "refined_bins": os.path.join(output_dir, "refined_bins"), # Maybe keep for copying final refined versions? Or remove. Let's remove for now.
        "checkm_results": os.path.join(output_dir, "checkm_results"),
        "figures": os.path.join(output_dir, "figures"),
        "final_bins": os.path.join(output_dir, "final_bins"),
        "temp": os.path.join(output_dir, "temp") # Keep for intermediate binning_refiner outputs and inputs
    }

    logger.info("Setting up output directories...")
    for name, path in subdirs.items():
         # Only create if not skipping relevant step or if it's essential (like temp)
        should_create = True
        if name == "checkm_results" and args.skip_checkm: should_create = False
        #if name == "refined_bins" and args.skip_refinement: should_create = False # Removed this subdir
        if name == "final_bins" and args.skip_consolidation: should_create = False
        if name == "figures" and args.skip_visualization: should_create = False
        # Always create temp dir unless keep_temp is false AND we are skipping refinement (no temp needed)
        if name == "temp" and not args.keep_temp and args.skip_refinement: should_create = False


        if should_create:
             logger.debug(f"Creating directory: {path}")
             os.makedirs(path, exist_ok=True)
        else:
            logger.debug(f"Skipping creation of directory: {path}")


    return subdirs


def clean_bin_folder_name(folder_path, bin_folder_to_name=None):
    """Clean bin folder names to create better combination/binner names."""
    if bin_folder_to_name and folder_path in bin_folder_to_name:
        return bin_folder_to_name[folder_path]

    # Extract basename
    name = os.path.basename(folder_path)

    # Remove common suffixes like '_bins', '.bins', '/bins'
    name = re.sub(r'[_\./]?bin[s]?$', '', name, flags=re.IGNORECASE)

    # Remove any other non-alphanumeric characters except underscores and hyphens
    name = re.sub(r'[^a-zA-Z0-9_\-]', '_', name)

    # Avoid names starting with numbers or hyphens if they might cause issues
    if name and (name[0].isdigit() or name[0] == '-'):
        name = f"b_{name}"

    # Ensure name is not empty
    if not name:
        name = f"binner_{hashlib.md5(folder_path.encode()).hexdigest()[:6]}" # Fallback to hash

    return name

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MODIFIED FUNCTION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def generate_bin_combinations(bin_folders, bin_folder_to_name=None, max_combinations=None):
    """Generate combinations of bin folders (size >= 2) up to a maximum number."""
    n = len(bin_folders)
    all_combinations = []
    min_size = 2 # Hardcoded minimum size for combinations

    if n < min_size:
         logger.warning(f"Cannot generate combinations of size >= {min_size} with only {n} input folders.")
         return [], {} # Return empty list and mapping


    # Start with largest combinations and work down to min_size
    sizes_to_generate = list(range(n, min_size - 1, -1))

    logger.debug(f"Generating combinations for sizes: {sizes_to_generate}")

    # Create sanitized names for bin folders
    folder_to_sanitized_name = {folder: clean_bin_folder_name(folder, bin_folder_to_name) for folder in bin_folders}
    sanitized_names_list = [folder_to_sanitized_name[folder] for folder in bin_folders] # Maintain order

    for size in sizes_to_generate:
        if size == 0: continue # Skip size 0 (shouldn't happen with range)
        size_combos = []
        for combo_indices in itertools.combinations(range(len(bin_folders)), size):
            combo_paths = [bin_folders[i] for i in combo_indices]
            combo_sanitized_names = [sanitized_names_list[i] for i in combo_indices]

            # Sort sanitized names within the combo for consistent naming
            combo_sanitized_names.sort()

            # Handle potential empty sanitized names (shouldn't happen with fallback)
            valid_sanitized_names = [n for n in combo_sanitized_names if n]
            if not valid_sanitized_names:
                 logger.warning(f"Skipping combination for indices {combo_indices} due to empty sanitized names.")
                 continue

            combo_info = {
                'size': size,
                'paths': combo_paths,
                'names': valid_sanitized_names, # Store the list of sanitized names
                'name': f"combo_size{size}_" + "_".join(valid_sanitized_names)
            }
            size_combos.append(combo_info)

        # Sort combinations within each size by names for reproducibility
        size_combos.sort(key=lambda x: x['name'])
        all_combinations.extend(size_combos)

    # Limit combinations if specified
    if max_combinations is not None and len(all_combinations) > max_combinations:
        logger.warning(f"Limiting to {max_combinations} combinations out of {len(all_combinations)} possible")
        all_combinations = all_combinations[:max_combinations]

    return all_combinations, folder_to_sanitized_name # Return mapping as well
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END MODIFIED FUNCTION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def force_create_folder(folder_to_create):
    """Creates a folder, removing the old one if it exists."""
    if os.path.isdir(folder_to_create):
        logger.debug(f"Removing existing directory: {folder_to_create}")
        shutil.rmtree(folder_to_create, ignore_errors=True)
    elif os.path.exists(folder_to_create): # Handle case where it exists but is not a dir
         logger.debug(f"Removing existing file/link: {folder_to_create}")
         try:
            os.remove(folder_to_create)
         except OSError as e:
            logger.warning(f"Could not remove existing file/link {folder_to_create}: {e}")


    logger.debug(f"Creating directory: {folder_to_create}")
    os.makedirs(folder_to_create)


def prepare_combination_input_dir(combination, output_dir, folder_to_sanitized_name):
    """Create the input directory structure expected by binning_refiner.py."""
    combo_input_base_dir = os.path.join(output_dir, combination['name'] + "_input") # Temp dir for this combo's input
    force_create_folder(combo_input_base_dir) # Ensure clean start

    for bin_folder_path in combination['paths']:
        binner_name = folder_to_sanitized_name.get(bin_folder_path)
        if not binner_name:
            logger.error(f"Could not find sanitized name for path: {bin_folder_path} in combination {combination['name']}")
            continue # Skip this folder if name mapping failed

        binner_input_dir = os.path.join(combo_input_base_dir, binner_name)
        os.makedirs(binner_input_dir, exist_ok=True)

        # Link or copy bin files into the binner-specific subdirectory
        files_linked_or_copied = 0
        try:
            for item in os.listdir(bin_folder_path):
                source_item_path = os.path.join(bin_folder_path, item)
                target_item_path = os.path.join(binner_input_dir, item)
                is_fasta, _ = is_fasta_file(item)

                if os.path.isfile(source_item_path) and is_fasta:
                    if not os.path.exists(target_item_path):
                        try:
                            # Use relative symlink if possible, might be more portable?
                            # source_rel = os.path.relpath(source_item_path, start=binner_input_dir)
                            # os.symlink(source_rel, target_item_path)
                            # Sticking to absolute for now, less likely to break
                            os.symlink(source_item_path, target_item_path)
                            files_linked_or_copied += 1
                        except OSError as e:
                            logger.warning(f"Symlink failed for {item} in {binner_name} (Error: {e}), copying instead.")
                            try:
                                shutil.copy2(source_item_path, target_item_path)
                                files_linked_or_copied += 1
                            except Exception as copy_e:
                                logger.error(f"Copy also failed for {item} in {binner_name}: {copy_e}")

            logger.debug(f"Prepared input for {binner_name} in {combination['name']} with {files_linked_or_copied} files.")
            if files_linked_or_copied == 0:
                 logger.warning(f"No bin files were linked or copied for binner {binner_name} from {bin_folder_path}")

        except Exception as e:
             logger.error(f"Error processing bin folder {bin_folder_path} for combo {combination['name']}: {e}")
             # Decide how to handle: maybe skip this combo or raise error?
             # For now, just log and continue, but the combo might fail later.
             continue

    return combo_input_base_dir


def process_bin_combination(combination, base_output_dir, refiner_script, min_bin_size, folder_to_sanitized_name, keep_temp=False):
    """Process a single bin combination using binning_refiner."""
    combo_name = combination['name']
    # Use the main 'temp' directory within the overall output directory
    temp_base_dir = os.path.join(base_output_dir, "temp")
    os.makedirs(temp_base_dir, exist_ok=True)

    combo_input_dir = None # Initialize
    refined_output_dir = None # Initialize
    output_prefix = None # Initialize

    try:
        # 1. Prepare the input directory for binning_refiner.py
        combo_input_dir = prepare_combination_input_dir(combination, temp_base_dir, folder_to_sanitized_name)

        # 2. Define output path for binning_refiner.py (within temp)
        output_prefix = os.path.join(temp_base_dir, combo_name + "_refined_output") # Place binning_refiner output here

        # 3. Run binning_refiner
        logger.info(f"Running binning_refiner on {combo_name}")
        logger.debug(f"Input dir for binning_refiner: {combo_input_dir}")
        logger.debug(f"Output prefix for binning_refiner: {output_prefix}")

        command = [
            sys.executable,  # Use the current Python interpreter
            refiner_script,  # Direct path to binning_refiner.py
            "-i", combo_input_dir,
            "-o", output_prefix, # Output within temp dir
            "-p", combo_name, # Use combo name as prefix for refined bins
            "-m", str(min_bin_size),
            "-q", "True"  # Silence binning_refiner output (use "False" for debugging it)
        ]
        # Check if binning_refiner supports -gz or similar for compressed output?
        # Let's assume based on previous structure output it *might* produce .gz
        # We don't pass a flag here, just handle whatever it produces.

        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True, # Capture output as text
            check=False # Don't raise error immediately, check manually
        )

        # Log stdout/stderr for debugging if needed
        if result.returncode != 0 or result.stderr: # Check stderr even if return code is 0
             logger.error(f"Binning_refiner potentially failed for {combo_name}. Return Code: {result.returncode}")
             # Only log stderr/stdout if they contain something
             if result.stderr:
                  logger.error(f"Stderr:\n{result.stderr}")
             if result.stdout:
                  logger.error(f"Stdout (potential errors):\n{result.stdout}") # Log stdout too, maybe error is there
             # Continue to check for output dir, maybe it partially succeeded?
        else:
             # Log stdout only if debug is enabled and it's not empty
             if logger.isEnabledFor(logging.DEBUG) and result.stdout:
                  logger.debug(f"Binning_refiner stdout for {combo_name}:\n{result.stdout}")


        # 4. Check if the standard 'Refined' folder was created by binning_refiner
        refined_output_dir = os.path.join(output_prefix, "Refined") # Expected standard output folder name

        if not os.path.isdir(refined_output_dir):
            logger.error(f"Binning_refiner did not create expected output directory: {refined_output_dir}")
            # Log directory listing for debugging
            if os.path.isdir(output_prefix):
                 logger.error(f"Contents of output prefix dir '{output_prefix}': {os.listdir(output_prefix)}")
            else:
                 logger.error(f"Output prefix directory '{output_prefix}' does not exist.")
            # Cleanup input dir if failed
            if combo_input_dir and os.path.exists(combo_input_dir) and not keep_temp:
                 shutil.rmtree(combo_input_dir)
            # Also cleanup the potentially failed output prefix dir
            if output_prefix and os.path.exists(output_prefix) and not keep_temp:
                 shutil.rmtree(output_prefix)

            return None # Indicate failure


        # 5. Check if any bins were produced inside the 'Refined' directory
        bin_files = []
        try:
            for filename in os.listdir(refined_output_dir):
                filepath = os.path.join(refined_output_dir, filename)
                if os.path.isfile(filepath):
                    is_fasta, _ = is_fasta_file(filename)
                    if is_fasta:
                        bin_files.append(filepath) # Store full path
        except OSError as e:
             logger.error(f"Cannot list files in refined output directory {refined_output_dir}: {e}")
             # Cleanup and indicate failure
             if combo_input_dir and os.path.exists(combo_input_dir) and not keep_temp: shutil.rmtree(combo_input_dir)
             if output_prefix and os.path.exists(output_prefix) and not keep_temp: shutil.rmtree(output_prefix)
             return None


        if not bin_files:
            logger.warning(f"No refined bins (fasta files) found in the output directory {refined_output_dir} for {combo_name}")
            # Cleanup - keep output_prefix dir if keep_temp, else remove it and input dir
            if not keep_temp:
                if output_prefix and os.path.exists(output_prefix): shutil.rmtree(output_prefix)
                if combo_input_dir and os.path.exists(combo_input_dir): shutil.rmtree(combo_input_dir)
            return None # Indicate no bins


        # 6. Check for filename length (binning_refiner already uses prefix)
        output_bin_count = 0
        final_bin_paths = []
        for bin_filepath in bin_files:
             bin_filename = os.path.basename(bin_filepath)
             target_filename = bin_filename # Keep original name from binning_refiner by default

             # Check filename length before potential final copy/move
             if len(target_filename) > FILENAME_MAX_LEN:
                  # Generate a shorter name if too long
                  fname_no_ext, ext = os.path.splitext(bin_filename)
                  # Handle double extensions like .fasta.gz
                  if ext.lower() == '.gz':
                      inner_fname_no_ext, inner_ext = os.path.splitext(fname_no_ext)
                      if inner_ext.lower() in ['.fa', '.fasta', '.fna']:
                           ext = inner_ext + ext
                           fname_no_ext = inner_fname_no_ext

                  hashed_name = hashlib.md5(fname_no_ext.encode()).hexdigest()[:16]
                  new_target_filename = f"{combo_name}_hash_{hashed_name}{ext}"
                  logger.warning(f"Original refined bin name '{bin_filename}' for combo '{combo_name}' is too long. Renaming to '{new_target_filename}'.")

                  target_path_new = os.path.join(refined_output_dir, new_target_filename)
                  try:
                       os.rename(bin_filepath, target_path_new)
                       final_bin_paths.append(target_path_new) # Add new path
                       output_bin_count += 1
                  except OSError as e:
                       logger.error(f"Failed to rename long filename {bin_filename} to {new_target_filename}: {e}")
                       # Keep original path if rename fails? Or skip? Skipping for now.
                       continue

             else:
                  final_bin_paths.append(bin_filepath) # Add original path
                  output_bin_count += 1


        logger.info(f"Successfully processed {combo_name}, found {output_bin_count} refined bins in {refined_output_dir}")
        # Don't cleanup input dir yet if keep_temp is True
        if not keep_temp and combo_input_dir and os.path.exists(combo_input_dir):
             shutil.rmtree(combo_input_dir)

        # Return the path to the directory containing the actual refined bins
        return refined_output_dir

    except subprocess.CalledProcessError as e:
        # This block might not be reached if check=False, but keep for safety
        logger.error(f"Error running binning_refiner on {combo_name}: {e.stderr}")
        if combo_input_dir and os.path.exists(combo_input_dir) and not keep_temp: shutil.rmtree(combo_input_dir)
        if output_prefix and os.path.exists(output_prefix) and not keep_temp: shutil.rmtree(output_prefix)
        return None
    except Exception as e:
        logger.exception(f"Unexpected error processing combination {combo_name}: {str(e)}") # Use logger.exception for traceback
        if combo_input_dir and os.path.exists(combo_input_dir) and not keep_temp: shutil.rmtree(combo_input_dir)
        if output_prefix and os.path.exists(output_prefix) and not keep_temp: shutil.rmtree(output_prefix)
        return None


def process_combinations_in_parallel(combinations, base_output_dir, refiner_script, min_bin_size, threads, batch_size, folder_to_sanitized_name, keep_temp):
    """Process bin combinations in parallel batches."""
    all_refined_output_dirs = [] # Store paths to the *output* directories containing refined bins
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
                    base_output_dir, # Pass main output dir
                    refiner_script,
                    min_bin_size,
                    folder_to_sanitized_name,
                    keep_temp # Pass keep_temp flag
                ): combo for combo in batch
            }


            for future in concurrent.futures.as_completed(future_to_combo):
                combo = future_to_combo[future]
                try:
                    # future.result() returns the path to the 'Refined' dir or None
                    refined_dir = future.result()
                    if refined_dir:
                        all_refined_output_dirs.append(refined_dir)
                    else:
                         # Log failure here as well, as process_bin_combination might return None without exception
                         logger.warning(f"Combination {combo['name']} did not produce a valid output directory.")
                except Exception as e:
                    # Exception was logged in process_bin_combination, maybe add combo name here
                    logger.error(f"Failure recorded for combination {combo['name']} during parallel processing: {e}")

    return all_refined_output_dirs # Return list of paths to refined output dirs


def collect_all_bins(bin_dirs_or_files):
    """Collect all bin files from multiple directories or directly from a list of files."""
    all_bin_paths = []
    logger.debug(f"Collecting bins from sources: {bin_dirs_or_files}")

    for item in bin_dirs_or_files:
        item_path = os.path.abspath(os.path.expanduser(item))
        if os.path.isdir(item_path):
            logger.debug(f"Scanning directory: {item_path}")
            try:
                for filename in os.listdir(item_path):
                    filepath = os.path.join(item_path, filename)
                    if os.path.isfile(filepath):
                         is_fasta, _ = is_fasta_file(filename)
                         if is_fasta:
                             all_bin_paths.append(filepath)
                             logger.debug(f"  Found bin file: {filepath}")
            except OSError as e:
                 logger.warning(f"Could not read directory {item_path}: {e}")
        elif os.path.isfile(item_path):
             is_fasta, _ = is_fasta_file(item)
             if is_fasta:
                  all_bin_paths.append(item_path)
                  logger.debug(f"Included bin file directly: {item_path}")
        else:
             logger.warning(f"Item is not a valid file or directory, skipping: {item}")

    # Ensure absolute paths and remove duplicates
    unique_absolute_paths = sorted(list(set(os.path.abspath(p) for p in all_bin_paths)))
    logger.info(f"Collected {len(unique_absolute_paths)} unique bin files.")
    return unique_absolute_paths


def run_unified_checkm2(all_bin_sources, checkm_output_dir, db_path, threads, args_output_dir, folder_to_sanitized_name):
    """Run CheckM2 once on all collected bin files (original and refined)."""
    checkm_results_df = pd.DataFrame() # Initialize empty DataFrame
    bin_mapping = {} # Will map input bin path to its source info

    # 1. Collect all bin file paths (absolute paths)
    # all_bin_sources will be a list containing original bin folders and refined bin output folders
    all_bin_filepaths = collect_all_bins(all_bin_sources)

    if not all_bin_filepaths:
        logger.warning("No bin files found to run CheckM2 on.")
        return checkm_results_df, bin_mapping # Return empty results


    # 2. Prepare bin mapping (maps absolute path to source set name and original filename)
    # Create reverse mapping for easier lookup
    sanitized_name_to_folder = {v: k for k, v in folder_to_sanitized_name.items()}

    for bin_path in all_bin_filepaths:
        bin_id = os.path.basename(bin_path)
        bin_set_name = "Unknown" # Default
        parent_dir = os.path.dirname(bin_path)
        parent_parent_dir = os.path.dirname(parent_dir)

        # Strategy:
        # a. Check if parent dir matches an original input folder path.
        # b. Check if parent dir is 'Refined' and parent-parent contains '_refined_output'.
        # c. Fallback to cleaned name of parent dir.

        found_match = False
        # a. Check against original input folders
        for orig_folder, san_name in folder_to_sanitized_name.items():
             if parent_dir == orig_folder:
                 bin_set_name = san_name
                 found_match = True
                 break

        # b. Check if it's a refined bin
        if not found_match and os.path.basename(parent_dir) == "Refined":
             parent_parent_base = os.path.basename(parent_parent_dir)
             if parent_parent_base.endswith("_refined_output"):
                  combo_name = parent_parent_base.replace("_refined_output", "")
                  bin_set_name = f"Refined_{combo_name}"
                  found_match = True

        # c. Fallback: Clean the parent dir name
        if not found_match:
             cleaned_parent_name = clean_bin_folder_name(parent_dir)
             # Check if this cleaned name corresponds to a known original binner
             if cleaned_parent_name in sanitized_name_to_folder:
                 bin_set_name = cleaned_parent_name
             else: # Truly unknown or intermediate dir? Use the cleaned name anyway.
                 bin_set_name = cleaned_parent_name if cleaned_parent_name else os.path.basename(parent_dir)

        bin_mapping[bin_path] = {
             "bin_set": bin_set_name,
             "bin_id": bin_id, # Keep original filename as ID within set
             "original_path": bin_path # Store the path used by checkm2
         }
        logger.debug(f"Mapped bin path {bin_path} to bin_set '{bin_set_name}'")


    # 3. Run CheckM2 using the list of file paths
    os.makedirs(checkm_output_dir, exist_ok=True)
    checkm_cmd = [
        "checkm2", "predict",
        "--threads", str(threads),
        "--input" # Use the list of files directly
    ]
    checkm_cmd.extend(all_bin_filepaths) # Add all file paths
    checkm_cmd.extend([
        "--output-directory", checkm_output_dir,
        "--force", # Overwrite previous results if necessary
        "--database_path", db_path
        # No --extension needed when providing a list of files
    ])

    logger.info(f"Running unified CheckM2 on {len(all_bin_filepaths)} bins.")
    logger.debug(f"CheckM2 command (list input):\ncheckm2 predict --threads {threads} --input [list_of_{len(all_bin_filepaths)}_files] --output-directory {checkm_output_dir} --force --database_path {db_path}")

    try:
        result = subprocess.run(checkm_cmd, check=True, capture_output=True, text=True)
        logger.info("CheckM2 completed successfully.")
        logger.debug(f"CheckM2 stdout:\n{result.stdout}")
        if result.stderr: # Log stderr even on success, as it might contain warnings
            logger.info(f"CheckM2 stderr:\n{result.stderr}")

        # 4. Parse results
        quality_report_path = os.path.join(checkm_output_dir, "quality_report.tsv")
        if os.path.exists(quality_report_path):
            checkm_results_df = pd.read_csv(quality_report_path, sep='\t')
            # Rename 'Name' column to 'original_path' to match our mapping key
            if 'Name' in checkm_results_df.columns:
                 checkm_results_df = checkm_results_df.rename(columns={'Name': 'original_path'})
            else:
                 logger.warning("CheckM2 quality report missing 'Name' column. Cannot map results reliably.")
                 # Attempt to proceed without mapping if column missing
                 return checkm_results_df, bin_mapping


            logger.info(f"Parsed CheckM2 results from {quality_report_path}")
            # Add bin_set and bin_id from mapping
            mapping_df = pd.DataFrame.from_dict(bin_mapping, orient='index').reset_index(drop=True)
            # Merge requires the common column 'original_path' to exist and match
            if 'original_path' in checkm_results_df.columns:
                 # Ensure paths match exactly (e.g., absolute vs relative)
                 # CheckM2 output 'Name' (now 'original_path') should be the absolute path we passed in
                 checkm_results_df = pd.merge(checkm_results_df, mapping_df[['original_path', 'bin_set', 'bin_id']], on='original_path', how='left')
                 # Check for merge failures (NaNs in added columns)
                 if checkm_results_df['bin_set'].isnull().any():
                      num_unmapped = checkm_results_df['bin_set'].isnull().sum()
                      logger.warning(f"{num_unmapped} CheckM2 results could not be mapped back to bin sets. Check paths in 'quality_report.tsv' and 'bin_mapping_path_to_info.pkl'.")
                      logger.debug(f"Example unmapped path from CheckM2: {checkm_results_df[checkm_results_df['bin_set'].isnull()]['original_path'].iloc[0]}")

            else:
                 logger.warning("Could not find 'original_path' column in CheckM2 output after renaming 'Name'. Cannot merge mapping details.")

        else:
            logger.error(f"CheckM2 quality report not found at: {quality_report_path}")

        # Save bin mapping for reference
        mapping_file = os.path.join(checkm_output_dir, "bin_mapping_path_to_info.pkl")
        try:
             with open(mapping_file, "wb") as f:
                 pickle.dump(bin_mapping, f)
             logger.info(f"Saved bin path mapping to {mapping_file}")
        except Exception as e:
             logger.error(f"Could not save bin mapping: {e}")


    except subprocess.CalledProcessError as e:
        logger.error(f"Error running CheckM2: {e}")
        logger.error(f"Stderr:\n{e.stderr}")
        logger.error(f"Stdout:\n{e.stdout}")
        # Return empty DataFrame if CheckM2 fails
        return pd.DataFrame(), bin_mapping
    except Exception as e:
         logger.exception(f"An unexpected error occurred during CheckM2 run: {e}")
         return pd.DataFrame(), bin_mapping


    return checkm_results_df, bin_mapping # Return dataframe and mapping


def get_contig_to_bin_mappings(bin_dirs_or_files, folder_to_sanitized_name):
    """Read all bins and return contig mappings and lengths."""
    contig_to_bins = defaultdict(list)
    contig_lengths = {}
    bin_info = {} # Store info about each bin (path, set_name, length, unique_bin_id)
    bin_sets_found = set()

    all_bin_filepaths = collect_all_bins(bin_dirs_or_files)

    logger.info(f"Processing {len(all_bin_filepaths)} bins for contig mapping...")
    sanitized_name_to_folder = {v: k for k, v in folder_to_sanitized_name.items()}


    for bin_path in all_bin_filepaths:
        try:
            # Attempt to determine the bin set name robustly
            binner_name = "Unknown"
            parent_dir = os.path.dirname(bin_path)
            parent_parent_dir = os.path.dirname(parent_dir) # Directory above the bin's directory

            # Strategy (same as in run_unified_checkm2 mapping):
            # 1. Check if parent dir path matches an original input folder path.
            # 2. Check if parent dir is 'Refined' and infer from parent-parent.
            # 3. Use cleaned name of parent dir as fallback.

            found_match = False
            parent_dir_base = os.path.basename(parent_dir)

            # 1. Check if the parent directory path matches an original input folder path
            for orig_folder, san_name in folder_to_sanitized_name.items():
                 if parent_dir == orig_folder:
                     binner_name = san_name
                     found_match = True
                     break

            # 2. If not matched to an original input, check if it's a refined bin
            if not found_match and parent_dir_base == "Refined":
                 parent_parent_base = os.path.basename(parent_parent_dir)
                 if parent_parent_base.endswith("_refined_output"):
                      combo_name = parent_parent_base.replace("_refined_output", "")
                      binner_name = f"Refined_{combo_name}"
                      found_match = True

            # 3. Fallback: Use the cleaned name of the immediate parent directory
            if not found_match:
                 cleaned_parent_name = clean_bin_folder_name(parent_dir)
                 # Check if this cleaned name corresponds to a known original binner
                 if cleaned_parent_name in sanitized_name_to_folder:
                     binner_name = cleaned_parent_name
                     found_match = True # Mark as found even if fallback
                 else: # Truly unknown or intermediate dir? Use the cleaned name anyway.
                     binner_name = cleaned_parent_name if cleaned_parent_name else parent_dir_base # Use original base name if cleaning fails


            bin_sets_found.add(binner_name)
            # Use a unique ID combining the inferred set name and the filename
            # Ensure filename itself is somewhat sanitized for use in ID
            safe_filename = re.sub(r'[^a-zA-Z0-9_\-\.]', '_', os.path.basename(bin_path))
            bin_id = f"{binner_name}:{safe_filename}"
            total_length = 0

            with open_fasta(bin_path) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    contig_id = str(record.id)
                    contig_len = len(record.seq)
                    contig_to_bins[contig_id].append(bin_id)
                    # Store length only once, warn if inconsistent?
                    if contig_id in contig_lengths and contig_lengths[contig_id] != contig_len:
                         logger.warning(f"Contig {contig_id} has inconsistent lengths ({contig_lengths[contig_id]} vs {contig_len}) in different bins. Using first length found.")
                    elif contig_id not in contig_lengths:
                        contig_lengths[contig_id] = contig_len

                    total_length += contig_len
            bin_info[bin_id] = {"path": bin_path, "length": total_length, "bin_set": binner_name}
            logger.debug(f"Processed bin {bin_id} from set {binner_name}, length {total_length}")

        except FileNotFoundError:
             logger.error(f"Bin file not found during contig mapping: {bin_path}")
        except Exception as e:
            logger.error(f"Error reading bin {bin_path}: {e}")

    logger.info(f"Found {len(contig_lengths)} unique contigs across {len(bin_info)} bins.")
    logger.info(f"Found contigs from {len(bin_sets_found)} distinct bin sets: {sorted(list(bin_sets_found))}")
    return contig_to_bins, contig_lengths, bin_info, list(bin_sets_found)


def calculate_bin_similarity(contig_to_bins, contig_lengths, bin_ids):
    """Calculate Jaccard similarity between pairs of bins based on shared contigs."""
    num_bins = len(bin_ids)
    if num_bins == 0:
        logger.warning("No bins provided for similarity calculation.")
        return np.array([]).reshape(0,0), np.array([]).reshape(0,0)

    similarity_matrix = np.zeros((num_bins, num_bins))
    bin_contigs = defaultdict(set)

    # Create a mapping from bin_id to its set of contigs
    logger.debug(f"Building contig sets for {num_bins} bins...")
    contigs_in_scope = set(contig_lengths.keys()) # Only consider contigs we have length info for
    for contig, bins in contig_to_bins.items():
        if contig not in contigs_in_scope: continue # Skip contigs without length info
        for bin_id in bins:
            # Only consider bins that are in the provided bin_ids list
            if bin_id in bin_ids:
                bin_contigs[bin_id].add(contig)

    # Calculate pairwise Jaccard similarity
    logger.debug(f"Calculating pairwise similarity matrix ({num_bins}x{num_bins})...")
    bin_id_map = {bin_id: i for i, bin_id in enumerate(bin_ids)} # Faster lookup

    for i in range(num_bins):
        bin1_id = bin_ids[i]
        contigs1 = bin_contigs.get(bin1_id, set())
        len1 = sum(contig_lengths.get(c, 0) for c in contigs1)

        for j in range(i, num_bins):
            bin2_id = bin_ids[j]
            contigs2 = bin_contigs.get(bin2_id, set())
            len2 = sum(contig_lengths.get(c, 0) for c in contigs2)


            intersection_contigs = contigs1.intersection(contigs2)
            union_contigs = contigs1.union(contigs2)

            # Weighted Jaccard by contig length
            intersection_len = sum(contig_lengths.get(c, 0) for c in intersection_contigs)
            # Calculate union length carefully to avoid double counting intersection
            # union_len = len1 + len2 - intersection_len
            # Alternative: sum lengths of union set
            union_len = sum(contig_lengths.get(c, 0) for c in union_contigs)


            if union_len > 0:
                similarity = intersection_len / union_len
            else:
                # Handle case where both bins are empty or only contain contigs not in contig_lengths
                similarity = 1.0 if i == j else 0.0 # Self similarity is 1, else 0

            idx1 = bin_id_map[bin1_id]
            idx2 = bin_id_map[bin2_id]
            similarity_matrix[idx1, idx2] = similarity
            similarity_matrix[idx2, idx1] = similarity # Ensure symmetry

    # Convert similarity to distance for clustering
    distance_matrix = 1 - similarity_matrix
    # Ensure distance is non-negative and handle potential floating point issues
    distance_matrix[distance_matrix < 0] = 0
    np.fill_diagonal(distance_matrix, 0) # Ensure diagonal is exactly zero

    logger.debug("Similarity calculation finished.")
    return distance_matrix, similarity_matrix


def cluster_bins(distance_matrix, bin_ids, min_cluster_size=2, min_samples=1):
    """Cluster bins using HDBSCAN."""
    # Ensure the distance matrix is in condensed format if necessary, or squareform
    # HDBSCAN takes a square distance matrix with metric='precomputed'
    num_bins = len(bin_ids)
    if distance_matrix.shape != (num_bins, num_bins):
         logger.error(f"Distance matrix shape {distance_matrix.shape} mismatch with number of bins {num_bins}.")
         return None, None

    # Check distance matrix properties before clustering
    logger.debug(f"Distance matrix properties before HDBSCAN:")
    logger.debug(f"  Shape: {distance_matrix.shape}")
    logger.debug(f"  Contains NaN: {np.isnan(distance_matrix).any()}")
    logger.debug(f"  Contains Inf: {np.isinf(distance_matrix).any()}")
    if distance_matrix.size > 0: # Avoid error on empty matrix
        logger.debug(f"  Min distance: {np.min(distance_matrix):.4f}, Max distance: {np.max(distance_matrix):.4f}")
        # Check symmetry
        is_symmetric = np.allclose(distance_matrix, distance_matrix.T)
        logger.debug(f"  Is symmetric: {is_symmetric}")
        if not is_symmetric:
             logger.warning("Distance matrix is not symmetric! HDBSCAN might behave unexpectedly.")
    else:
        logger.debug("  Distance matrix is empty.")
        return {}, [] # Return empty clusters and singletons if matrix is empty

    # Check for NaNs/Infs which cause HDBSCAN errors
    if np.isnan(distance_matrix).any() or np.isinf(distance_matrix).any():
        logger.error("Distance matrix contains NaN or Inf values. Cannot perform HDBSCAN clustering.")
        # Try to identify problematic bins?
        nan_rows, nan_cols = np.where(np.isnan(distance_matrix))
        inf_rows, inf_cols = np.where(np.isinf(distance_matrix))
        problem_indices = set(nan_rows) | set(nan_cols) | set(inf_rows) | set(inf_cols)
        problem_bins = [bin_ids[i] for i in problem_indices if i < len(bin_ids)]
        logger.error(f"Problematic bins (indices {list(problem_indices)}): {problem_bins}")
        return None, None

    # Explicitly check min_cluster_size before passing to HDBSCAN
    if min_cluster_size <= 1:
        logger.error(f"min_cluster_size ({min_cluster_size}) must be greater than 1 for HDBSCAN.")
        return None, None

    logger.debug(f"Running HDBSCAN with min_cluster_size={min_cluster_size}, min_samples={min_samples}, metric='precomputed'")

    try:
        # Ensure dtype is float64 for HDBSCAN compatibility
        distance_matrix_float = distance_matrix.astype(np.float64)

        clusterer = hdbscan.HDBSCAN(metric='precomputed',
                                    min_cluster_size=min_cluster_size,
                                    min_samples=min_samples,
                                    allow_single_cluster=True, # Allow single clusters
                                    gen_min_span_tree=True) # Helps with certain diagnostics if needed later
        cluster_labels = clusterer.fit_predict(distance_matrix_float) # Use float matrix
    except ValueError as e:
         # Catch potential ValueErrors from HDBSCAN (e.g., related to matrix properties)
         logger.error(f"HDBSCAN clustering failed with ValueError: {e}")
         return None, None
    except Exception as e:
         # Catch any other unexpected errors during clustering
         logger.exception(f"An unexpected error occurred during HDBSCAN: {e}")
         return None, None


    # Group bins by cluster label
    clusters = defaultdict(list)
    for bin_id, label in zip(bin_ids, cluster_labels):
        clusters[label].append(bin_id)

    # Separate singleton points (label -1)
    singleton_points = clusters.pop(-1, [])

    logger.info(f"HDBSCAN clustering resulted in {len(clusters)} clusters and {len(singleton_points)} singleton points.")
    if logger.isEnabledFor(logging.DEBUG):
         for label, bins in clusters.items():
              logger.debug(f"  Cluster {label}: {len(bins)} bins - {bins[:5]}...") # Log first few bins per cluster
         logger.debug(f"  Singleton points ({len(singleton_points)}): {singleton_points[:5]}...")

    return clusters, singleton_points


def validate_clusters(clusters, similarity_matrix, bin_ids_list, threshold=0.7):
    """Validate clusters: ensure all bins within a cluster meet the similarity threshold pairwise."""
    valid_clusters = {}
    split_bins = [] # Bins removed from invalid clusters
    new_singleton_clusters = {}
    singleton_counter = 0

    if not clusters:
         logger.info("No clusters found to validate.")
         return {}, []

    bin_id_to_index = {bin_id: idx for idx, bin_id in enumerate(bin_ids_list)}

    logger.info(f"Validating {len(clusters)} clusters with threshold {threshold}...")

    for label, bin_list in clusters.items():
        if not isinstance(bin_list, list): # Ensure it's a list
            logger.warning(f"Cluster {label} members is not a list: {type(bin_list)}. Skipping validation.")
            continue

        if len(bin_list) < 2:
            valid_clusters[label] = bin_list # Singletons are valid by definition
            continue

        # Create copy to modify if needed
        current_cluster_bins = list(bin_list)
        logger.debug(f"Validating cluster {label} (initial size {len(current_cluster_bins)})...")


        # Iteratively check and split until the cluster is stable
        iteration = 0
        while True:
            iteration += 1
            logger.debug(f"  Validation iteration {iteration} for cluster {label}")
            bins_to_remove_this_round = set()
            pairwise_checked = False # Flag to see if any pairs were checked

            # Check all pairs in the current list of bins for this cluster
            indices_in_cluster = [bin_id_to_index.get(b) for b in current_cluster_bins]
            # Filter out None indices if any bins weren't found in the master list
            valid_indices = [idx for idx in indices_in_cluster if idx is not None]
            if len(valid_indices) != len(current_cluster_bins):
                 missing_bins = [b for b, idx in zip(current_cluster_bins, indices_in_cluster) if idx is None]
                 logger.warning(f"Cluster {label}: Cannot find indices for bins: {missing_bins}. Proceeding with available bins.")
                 current_cluster_bins = [b for b, idx in zip(current_cluster_bins, indices_in_cluster) if idx is not None]
                 if len(current_cluster_bins) < 2: break # Stop if cluster becomes too small


            for i in range(len(current_cluster_bins)):
                for j in range(i + 1, len(current_cluster_bins)):
                    pairwise_checked = True
                    bin1_id = current_cluster_bins[i]
                    bin2_id = current_cluster_bins[j]
                    idx1 = bin_id_to_index.get(bin1_id) # Should exist now
                    idx2 = bin_id_to_index.get(bin2_id) # Should exist now

                    similarity = similarity_matrix[idx1, idx2]
                    if similarity < threshold:
                        logger.warning(f"  Cluster {label} invalid pair: {bin1_id} & {bin2_id} (Sim: {similarity:.3f} < {threshold})")
                        # Strategy: Remove the bin with lower average similarity to the rest of the current cluster members
                        sim1_avg = 0.0
                        sim2_avg = 0.0
                        count = 0
                        # Calculate average similarity excluding self and the other pair member
                        for k in range(len(current_cluster_bins)):
                             if k != i and k != j:
                                  bin_k_id = current_cluster_bins[k]
                                  idx_k = bin_id_to_index.get(bin_k_id)
                                  if idx_k is not None:
                                       sim1_avg += similarity_matrix[idx1, idx_k]
                                       sim2_avg += similarity_matrix[idx2, idx_k]
                                       count += 1
                        if count > 0:
                             sim1_avg /= count
                             sim2_avg /= count

                        # Decide which to remove
                        if sim1_avg < sim2_avg:
                             bins_to_remove_this_round.add(bin1_id)
                             logger.debug(f"    -> Marking {bin1_id} for removal (avg sim {sim1_avg:.3f} vs {sim2_avg:.3f})")
                        else: # Remove bin2 if scores are equal or bin2 is lower
                             bins_to_remove_this_round.add(bin2_id)
                             logger.debug(f"    -> Marking {bin2_id} for removal (avg sim {sim2_avg:.3f} vs {sim1_avg:.3f})")


            if bins_to_remove_this_round:
                 # Create new singletons for removed bins
                 for removed_bin in bins_to_remove_this_round:
                     singleton_counter += 1
                     new_label = f"split_{label}_{singleton_counter:03d}" # Use descriptive label
                     logger.info(f"Split bin {removed_bin} from cluster {label} into singleton cluster {new_label}")
                     new_singleton_clusters[new_label] = [removed_bin]
                     split_bins.append(removed_bin)

                 # Update current cluster list and re-check
                 current_cluster_bins = [b for b in current_cluster_bins if b not in bins_to_remove_this_round]
                 logger.debug(f"  Cluster {label} reduced to {len(current_cluster_bins)} bins. Revalidating.")
                 if len(current_cluster_bins) < 2:
                      logger.debug(f"  Cluster {label} now too small after iteration {iteration}.")
                      break # Stop if cluster becomes too small

            elif not pairwise_checked and len(current_cluster_bins) >=2 :
                 # This case might happen if indices were missing for all pairs
                 logger.warning(f"Could not perform pairwise validation for cluster {label} due to missing bin indices.")
                 # Treat as valid for now, but log warning.
                 break
            else:
                # No bins removed in this round, cluster is stable
                logger.debug(f"  Cluster {label} stable after iteration {iteration}.")
                break


        # After loop, decide what to do with the remaining bins
        if len(current_cluster_bins) >= 2:
            logger.debug(f"Validated cluster {label} with {len(current_cluster_bins)} bins.")
            valid_clusters[label] = current_cluster_bins
        elif len(current_cluster_bins) == 1:
             # Make the remaining bin a singleton
             singleton_counter += 1
             new_label = f"split_{label}_{singleton_counter:03d}"
             logger.info(f"Making remaining bin {current_cluster_bins[0]} from cluster {label} a singleton {new_label}")
             new_singleton_clusters[new_label] = current_cluster_bins
             split_bins.append(current_cluster_bins[0])
        else:
             logger.info(f"Cluster {label} dissolved completely after validation.")


    # Add the newly created singletons to the valid clusters
    valid_clusters.update(new_singleton_clusters)

    return valid_clusters, split_bins


def select_best_bin(cluster_bins, checkm_results_df, bin_info):
    """Select the best bin from a cluster based on CheckM quality."""
    best_bin = None
    best_score = -float('inf') # Initialize with negative infinity

    if not cluster_bins: # Handle empty cluster input
        logger.warning("Received empty cluster for bin selection.")
        return None

    # Use first bin ID for logging context if available
    cluster_id_str = f"Cluster starting with: {cluster_bins[0]}" if cluster_bins else "Empty Cluster"
    logger.debug(f"Selecting best bin for {cluster_id_str} (size: {len(cluster_bins)})")


    if checkm_results_df is None or checkm_results_df.empty:
        logger.warning(f"{cluster_id_str}: No CheckM results available. Falling back to largest bin.")
        # If no CheckM results, fall back to largest bin
        max_len = -1
        for bin_id in cluster_bins:
             info = bin_info.get(bin_id)
             if info and info.get('length', -1) > max_len:
                 max_len = info['length']
                 best_bin = bin_id
        if best_bin:
             logger.info(f"Selected largest bin by length: {best_bin} (Length: {max_len}) for {cluster_id_str}")
        else:
             logger.warning(f"Could not determine largest bin for {cluster_id_str}")
        return best_bin

    # Filter CheckM results for bins in the current cluster
    cluster_paths = [bin_info.get(b_id, {}).get('path') for b_id in cluster_bins if bin_info.get(b_id)]
    cluster_paths = [p for p in cluster_paths if p] # Remove None paths
    if not cluster_paths:
         logger.warning(f"Could not get paths for any bins in {cluster_id_str}. Falling back to largest bin.")
         # Fallback logic
         max_len = -1; best_bin = None # Reset
         for bin_id in cluster_bins:
              info = bin_info.get(bin_id)
              if info and info.get('length', -1) > max_len:
                  max_len = info['length']
                  best_bin = bin_id
         if best_bin: logger.info(f"Selected largest bin by length: {best_bin} (Length: {max_len}) for {cluster_id_str}")
         else: logger.warning(f"Could not determine largest bin for {cluster_id_str}")
         return best_bin


    cluster_checkm = checkm_results_df[checkm_results_df['original_path'].isin(cluster_paths)].copy()


    if cluster_checkm.empty:
         logger.warning(f"No CheckM results found for any bin in {cluster_id_str}. Falling back to largest bin.")
         # Fallback to largest bin
         max_len = -1; best_bin = None # Reset
         for bin_id in cluster_bins:
              info = bin_info.get(bin_id)
              if info and info.get('length', -1) > max_len:
                  max_len = info['length']
                  best_bin = bin_id
         if best_bin:
              logger.info(f"Selected largest bin by length: {best_bin} (Length: {max_len}) for {cluster_id_str}")
         else:
              logger.warning(f"Could not determine largest bin for {cluster_id_str}")
         return best_bin


    # Calculate score: Completion - 5 * Contamination
    if 'Completeness' in cluster_checkm.columns and 'Contamination' in cluster_checkm.columns:
        cluster_checkm['Completeness'] = pd.to_numeric(cluster_checkm['Completeness'], errors='coerce')
        cluster_checkm['Contamination'] = pd.to_numeric(cluster_checkm['Contamination'], errors='coerce')
        # Drop rows where conversion failed or quality is missing
        cluster_checkm = cluster_checkm.dropna(subset=['Completeness', 'Contamination'])

        if not cluster_checkm.empty:
             cluster_checkm['Score'] = cluster_checkm['Completeness'] - 5 * cluster_checkm['Contamination']
             # Add original length to break ties
             cluster_checkm['Length'] = cluster_checkm['original_path'].map(lambda p: bin_info.get(next((bid for bid, info in bin_info.items() if info.get('path') == p), None), {}).get('length', 0))

             # Handle ties: prefer lower contamination, then higher completeness, then greater length
             cluster_checkm = cluster_checkm.sort_values(
                  by=['Score', 'Contamination', 'Completeness', 'Length'],
                  ascending=[False, True, False, False] # Max Score, Min Contam, Max Comp, Max Length
             )
             best_row = cluster_checkm.iloc[0]
             best_bin_path = best_row['original_path']
             best_score = best_row['Score']

             # Find the bin_id corresponding to this path
             for b_id, info in bin_info.items():
                  if info.get('path') == best_bin_path:
                      best_bin = b_id
                      break # Found the match

             if best_bin:
                  logger.info(f"Best bin selected for {cluster_id_str}: {best_bin} (Score: {best_score:.2f}, C:{best_row['Completeness']:.2f}%, X:{best_row['Contamination']:.2f}%)")
                  return best_bin
             else:
                 logger.error(f"Could not map best bin path {best_bin_path} back to bin_id in {cluster_id_str}. Critical error.")
                 # Fallback just in case
                 max_len = -1; best_bin = None # Reset
                 for bin_id_fb in cluster_bins: info_fb = bin_info.get(bin_id_fb); # Shorten line
                 if info_fb and info_fb.get('length', -1) > max_len: max_len = info_fb['length']; best_bin = bin_id_fb
                 logger.warning(f"Falling back to largest bin by length for {cluster_id_str}: {best_bin}")
                 return best_bin

        else:
             logger.warning(f"Could not calculate scores (numeric conversion failed or no valid data) for bins in {cluster_id_str}. Falling back to largest bin.")
             # Fallback logic
             max_len = -1; best_bin = None # Reset
             for bin_id in cluster_bins: info = bin_info.get(bin_id)
             if info and info.get('length', -1) > max_len: max_len = info['length']; best_bin = bin_id
             if best_bin: logger.info(f"Selected largest bin by length: {best_bin} (Length: {max_len}) for {cluster_id_str}")
             else: logger.warning(f"Could not determine largest bin for {cluster_id_str}")
             return best_bin

    else:
         logger.warning(f"{cluster_id_str}: CheckM results missing 'Completeness' or 'Contamination' column. Falling back.")
         # Fallback logic
         max_len = -1; best_bin = None # Reset
         for bin_id in cluster_bins: info = bin_info.get(bin_id)
         if info and info.get('length', -1) > max_len: max_len = info['length']; best_bin = bin_id
         if best_bin: logger.info(f"Selected largest bin by length: {best_bin} (Length: {max_len}) for {cluster_id_str}")
         else: logger.warning(f"Could not determine largest bin for {cluster_id_str}")
         return best_bin


def resolve_ambiguous_contigs(best_bins, contig_to_bins, bin_info, keep_ambiguous=False, remove_ambiguous=False):
    """Handle contigs appearing in multiple selected best bins."""
    contigs_in_final_bins = defaultdict(list)
    final_bins_contigs = defaultdict(set)

    # Map contigs to the final selected bins
    for bin_id in best_bins:
        bin_path = bin_info.get(bin_id, {}).get('path')
        if not bin_path:
             logger.warning(f"Skipping bin {bin_id} in ambiguity check: path not found in bin_info.")
             continue
        try:
            with open_fasta(bin_path) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    contig_id = str(record.id)
                    contigs_in_final_bins[contig_id].append(bin_id)
                    final_bins_contigs[bin_id].add(contig_id)
        except FileNotFoundError:
             logger.error(f"Error reading final selected bin {bin_path}: File not found.")
        except Exception as e:
            logger.error(f"Error reading final selected bin {bin_path}: {e}")

    contigs_to_remove = defaultdict(set) # {bin_id: set_of_contigs_to_remove}
    # Find contigs present in more than one selected bin
    ambiguous_contigs = {c: list(set(bins)) for c, bins in contigs_in_final_bins.items() if len(set(bins)) > 1}


    if not ambiguous_contigs:
        logger.info("No ambiguous contigs found among selected best bins.")
        return final_bins_contigs # Return original contig sets

    logger.info(f"Found {len(ambiguous_contigs)} ambiguous contigs present in multiple selected bins.")
    logger.debug(f"Ambiguous contigs (first 5): {dict(itertools.islice(ambiguous_contigs.items(), 5))}")


    if remove_ambiguous:
        logger.info("Removing ambiguous contigs from all bins they appear in (--remove-ambiguous).")
        for contig_id, bins in ambiguous_contigs.items():
            for bin_id in bins:
                contigs_to_remove[bin_id].add(contig_id)
    elif not keep_ambiguous: # Default behavior if neither flag is set
         # Defaulting to *removing* ambiguous contigs if neither flag is set, for safety.
         logger.warning("Defaulting to removing ambiguous contigs as neither --keep-ambiguous nor --remove-ambiguous was specified.")
         for contig_id, bins in ambiguous_contigs.items():
            for bin_id in bins:
                contigs_to_remove[bin_id].add(contig_id)

    elif keep_ambiguous:
         logger.info("Keeping ambiguous contigs in all bins (--keep-ambiguous specified).")
         # Do nothing, contigs remain in all lists
         pass


    # Update final_bins_contigs by removing specified contigs
    if contigs_to_remove:
         bins_made_empty_count = 0
         updated_final_bins_contigs = {}
         contigs_removed_count = 0
         contigs_actually_removed = set()

         for bin_id, contigs_set in final_bins_contigs.items():
             original_size = len(contigs_set)
             to_remove_for_this_bin = contigs_to_remove.get(bin_id, set())
             # Find contigs that are actually present and marked for removal
             removed_in_this_bin = contigs_set.intersection(to_remove_for_this_bin)
             num_removed_this_bin = len(removed_in_this_bin)
             contigs_removed_count += num_removed_this_bin
             contigs_actually_removed.update(removed_in_this_bin)

             # Remove them
             contigs_set.difference_update(to_remove_for_this_bin)


             if contigs_set: # Only keep non-empty bins
                  updated_final_bins_contigs[bin_id] = contigs_set
             elif original_size > 0: # Only count if it was non-empty before
                 bins_made_empty_count += 1
                 logger.debug(f"Bin {bin_id} became empty after removing ambiguous contigs.")

         final_bins_contigs = updated_final_bins_contigs
         logger.info(f"Removed {contigs_removed_count} ambiguous contig occurrences ({len(contigs_actually_removed)} unique contigs) across bins.")
         if bins_made_empty_count > 0:
              logger.info(f"{bins_made_empty_count} bins became empty after removing ambiguous contigs.")


    return final_bins_contigs


def write_final_bins(final_bins_contigs, bin_info, output_dir):
    """Write the final set of bins after ambiguity resolution."""
    logger.info(f"Attempting to write {len(final_bins_contigs)} final bins to {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    written_count = 0

    # Keep track of original paths to avoid reading the same file multiple times
    contig_cache = {} # {original_path: {contig_id: record}}

    for bin_id, contigs_set in final_bins_contigs.items():
        if not contigs_set:
             logger.debug(f"Skipping bin_id {bin_id} as its contig set is empty.")
             continue # Skip empty bins

        # Sanitize bin_id for filename
        safe_bin_filename = re.sub(r'[^a-zA-Z0-9_\-\.]', '_', bin_id) # Allow dots and hyphens
        # Ensure it ends with .fasta (or keeps original extension if present)
        orig_bin_path = bin_info.get(bin_id, {}).get('path', '')
        if not orig_bin_path:
             logger.warning(f"Could not find original path for bin_id {bin_id}. Cannot determine output extension. Using .fasta.")
             orig_ext = ".fasta"
        else:
            _, orig_ext = os.path.splitext(orig_bin_path)
            # Handle double extensions like .fasta.gz
            if orig_ext.lower() == '.gz':
                 inner_fname_no_ext, inner_ext = os.path.splitext(os.path.splitext(orig_bin_path)[0])
                 if inner_ext.lower() in ['.fa', '.fasta', '.fna']:
                      orig_ext = inner_ext + orig_ext

        # Add extension if not present in sanitized name, otherwise keep existing
        has_fasta_ext = any(safe_bin_filename.lower().endswith(ext) for ext in ['.fa', '.fasta', '.fna'])
        has_gz_ext = safe_bin_filename.lower().endswith('.gz')

        if has_gz_ext and not has_fasta_ext: # e.g., mybin.gz -> mybin.fasta.gz
             safe_bin_filename = os.path.splitext(safe_bin_filename)[0] + ".fasta.gz"
        elif not has_fasta_ext: # Add original ext (could be .fasta.gz) or default to .fasta
             safe_bin_filename += orig_ext if orig_ext else ".fasta"


        final_bin_path = os.path.join(output_dir, safe_bin_filename)

        original_bin_path = bin_info.get(bin_id, {}).get('path')
        if not original_bin_path:
             logger.warning(f"Could not find original path for bin_id {bin_id} during writing. Skipping.")
             continue

        try:
            # Determine write mode based on final path extension
            write_mode = 'w'
            output_handle_opener = open # Default to standard open

            if final_bin_path.lower().endswith(".gz"):
                write_mode = 'wt' # Text mode for gzip
                output_handle_opener = gzip.open
                logger.debug(f"Writing gzipped output to: {final_bin_path}")
            else:
                 logger.debug(f"Writing plain text output to: {final_bin_path}")


            # Load original contigs if not cached
            if original_bin_path not in contig_cache:
                  logger.debug(f"Caching contigs from: {original_bin_path}")
                  contig_cache[original_bin_path] = {}
                  try:
                      # Use open_fasta to handle potential gzip input
                      with open_fasta(original_bin_path, 'r') as infile:
                           for record in SeqIO.parse(infile, "fasta"):
                               contig_cache[original_bin_path][str(record.id)] = record
                  except FileNotFoundError:
                       logger.error(f"Failed to cache - source file not found: {original_bin_path}")
                       continue # Skip this bin if source is missing
                  except Exception as parse_err:
                       logger.error(f"Failed to parse or cache {original_bin_path}: {parse_err}")
                       # Continue to next bin if caching fails for this source
                       continue

            # Write selected contigs using the correct handle
            with output_handle_opener(final_bin_path, write_mode, encoding='utf-8') as outfile:
                 contigs_written_count = 0
                 missing_contigs_in_cache = 0
                 for contig_id in sorted(list(contigs_set)): # Write in sorted order
                      record = contig_cache[original_bin_path].get(contig_id)
                      if record:
                           SeqIO.write(record, outfile, "fasta")
                           contigs_written_count += 1
                      else:
                           logger.warning(f"Contig {contig_id} expected in {bin_id} (source {original_bin_path}) but not found in cache.")
                           missing_contigs_in_cache +=1

                 if contigs_written_count > 0:
                      logger.debug(f"Wrote {contigs_written_count} contigs to {final_bin_path}" + (f" ({missing_contigs_in_cache} missing from cache)" if missing_contigs_in_cache else ""))
                      written_count += 1
                 else:
                      logger.warning(f"Final bin file {final_bin_path} was not written because no contigs were found/cached for it.")
                      # Attempt to remove the empty file if it was created (handle needs closing first)
                      # outfile handle is closed automatically by 'with' statement here.
                      try:
                           os.remove(final_bin_path)
                           logger.debug(f"Removed empty file: {final_bin_path}")
                      except Exception as remove_err:
                            logger.warning(f"Failed to remove empty file {final_bin_path}: {remove_err}")


        except Exception as e:
             logger.error(f"Error writing final bin {final_bin_path}: {e}")

    logger.info(f"Successfully wrote {written_count} final bins.")
    return written_count


def create_final_report(best_bins, checkm_results_df, bin_info, report_path):
     """Creates a TSV report summarizing the selected best bins and their quality."""
     report_data = []
     if not best_bins:
          logger.warning("No best bins selected, cannot create final report.")
          # Create empty report with headers
          pd.DataFrame(columns=[
                'Selected_BinID', 'Source_Set', 'Original_Path', 'Length',
                'Completeness', 'Contamination', 'Quality_Score'
            ]).to_csv(report_path, sep='\t', index=False)
          return


     for bin_id in best_bins:
          info = bin_info.get(bin_id, {})
          path = info.get('path', 'N/A')
          bin_set = info.get('bin_set', 'Unknown')
          length = info.get('length', 0)

          quality_info = {'Completeness': 'N/A', 'Contamination': 'N/A', 'Quality Score': 'N/A'}
          if checkm_results_df is not None and not checkm_results_df.empty and 'bin_id' in checkm_results_df.columns:
                # Match using the absolute path stored in bin_info
                checkm_row = checkm_results_df[checkm_results_df['bin_id'] == bin_id]
                if not checkm_row.empty:
                     row = checkm_row.iloc[0]
                     comp_val = row.get('Completeness', 'N/A')
                     cont_val = row.get('Contamination', 'N/A')
                     quality_info['Completeness'] = f"{comp_val:.2f}" if isinstance(comp_val, (int, float)) else comp_val
                     quality_info['Contamination'] = f"{cont_val:.2f}" if isinstance(cont_val, (int, float)) else cont_val

                     # Recalculate score if columns exist and are numeric
                     try:
                          comp = pd.to_numeric(comp_val, errors='coerce')
                          cont = pd.to_numeric(cont_val, errors='coerce')
                          if pd.notna(comp) and pd.notna(cont):
                               quality_info['Quality Score'] = f"{comp - 5 * cont:.2f}"
                          else:
                               # Handle cases where quality is 'N/A'
                               quality_info['Quality Score'] = 'N/A'
                     except Exception as e:
                          logger.debug(f"Could not calculate quality score for {bin_id}: {e}")
                          quality_info['Quality Score'] = 'N/A'


          report_data.append({
                'Selected_BinID': bin_id,
                'Source_Set': bin_set,
                'Original_Path': path,
                'Length': length,
                'Completeness': quality_info['Completeness'],
                'Contamination': quality_info['Contamination'],
                'Quality_Score': quality_info['Quality Score']
          })

     report_df = pd.DataFrame(report_data)
     # Sort report for consistency? Maybe by score or name?
     # report_df = report_df.sort_values(by='Selected_BinID')
     report_df.to_csv(report_path, sep='\t', index=False)
     logger.info(f"Final report saved to {report_path}")


def plot_bin_quality(best_bins, checkm_results_df, bin_info, plot_path):
     """Creates a scatter plot of final bin quality."""
     plot_data = []
     if checkm_results_df is None or checkm_results_df.empty or 'bin_id' not in checkm_results_df.columns:
         logger.warning("Cannot create quality plot: CheckM2 results are missing or incomplete.")
         return

     if not best_bins:
          logger.warning("No best bins selected, cannot create quality plot.")
          return


     for bin_id in best_bins:
         info = bin_info.get(bin_id, {})
         

         checkm_row = checkm_results_df[checkm_results_df['bin_id'] == bin_id]
         if not checkm_row.empty:
             row = checkm_row.iloc[0]
             try:
                 comp = pd.to_numeric(row.get('Completeness'), errors='coerce')
                 cont = pd.to_numeric(row.get('Contamination'), errors='coerce')
                 length = info.get('length', 0)
                 if pd.notna(comp) and pd.notna(cont):
                      plot_data.append({'Completeness': comp, 'Contamination': cont, 'Length': length, 'BinID': bin_id})
             except:
                  logger.warning(f"Could not parse quality for bin {bin_id} for plotting.")


     if not plot_data:
          logger.warning("No valid quality data to plot.")
          return

     plot_df = pd.DataFrame(plot_data)

     try:
          plt.figure(figsize=(10, 6))
          sns.scatterplot(data=plot_df, x='Completeness', y='Contamination', size='Length', hue='Length', sizes=(20, 200), legend="brief")
          plt.title('Final Bin Quality')
          plt.xlim(0, 105)
          plt.ylim(-1, max(15, plot_df['Contamination'].max() * 1.1 if not plot_df['Contamination'].empty else 15)) # Adjust ylim based on data
          plt.grid(True, linestyle='--', alpha=0.6)
          plt.savefig(plot_path, dpi=150, bbox_inches='tight')
          plt.close()
          logger.info(f"Bin quality plot saved to {plot_path}")
     except Exception as e:
          logger.error(f"Failed to create bin quality plot: {e}")
          # Attempt to close plot if it was opened
          try: plt.close()
          except: pass


def run_checkm2_on_final_bins(final_bins_dir, checkm_output_dir, db_path, threads):
     """Run CheckM2 specifically on the final consolidated bins."""
     logger.info(f"Running CheckM2 on final bins in directory: {final_bins_dir}")
     bin_files = collect_all_bins([final_bins_dir])
     if not bin_files:
         logger.info("No final bins found to run CheckM2 on.")
         return

     os.makedirs(checkm_output_dir, exist_ok=True)
     checkm_cmd = [
         "checkm2", "predict",
         "--threads", str(threads),
         "--input"
     ]
     checkm_cmd.extend(bin_files)
     checkm_cmd.extend([
         "--output-directory", checkm_output_dir,
         "--force",
         "--database_path", db_path
     ])

     logger.info(f"Running CheckM2 on {len(bin_files)} final bins.")
     logger.debug(f"Final CheckM2 command:\ncheckm2 predict --threads {threads} --input [list_of_{len(bin_files)}_files] --output-directory {checkm_output_dir} --force --database_path {db_path}")

     try:
         result = subprocess.run(checkm_cmd, check=True, capture_output=True, text=True)
         logger.info("CheckM2 on final bins completed successfully.")
         logger.debug(f"CheckM2 stdout:\n{result.stdout}")
         if result.stderr:
             logger.info(f"CheckM2 stderr:\n{result.stderr}")
         # Copy quality report to main directory for easier access?
         final_report_src = os.path.join(checkm_output_dir, "quality_report.tsv")
         final_report_dst = os.path.join(os.path.dirname(checkm_output_dir), "final_bins_checkm2_report.tsv")
         if os.path.exists(final_report_src):
             shutil.copy2(final_report_src, final_report_dst)
             logger.info(f"Copied final CheckM2 report to: {final_report_dst}")

     except subprocess.CalledProcessError as e:
         logger.error(f"Error running CheckM2 on final bins: {e}")
         logger.error(f"Stderr:\n{e.stderr}")
         logger.error(f"Stdout:\n{e.stdout}")
     except Exception as e:
         logger.exception(f"An unexpected error occurred during final CheckM2 run: {e}")


# ========================================
# Main Workflow
# ========================================
def main():
    args = parse_arguments()
    validate_inputs(args)
    output_dirs = setup_output_directory(args.output_dir, args) # Pass args here


    # --- Generate Binner Name Mapping ---
    binner_mapping = {} # folder_path -> binner_name (e.g., binner_001)
    if args.binner_names:
        provided_names = args.binner_names.split(',')
        for i, folder in enumerate(args.bin_folders):
             binner_mapping[folder] = provided_names[i]
    else:
        # Generate default names if not provided
        temp_name_counts = defaultdict(int)
        for i, folder in enumerate(args.bin_folders):
            base_name = clean_bin_folder_name(folder)
            # Use indexed name if cleaned name is generic, empty, or already used as a base
            if not base_name or base_name in ['bins', '_'] or temp_name_counts[base_name] > 0:
                 binner_name = f"binner_{i+1:03d}"
                 # Ensure this indexed name hasn't been manually provided and caused a clash
                 if any(name == binner_name for name in binner_mapping.values()):
                     # Extremely unlikely fallback, but possible
                     binner_name = f"binner_hash_{hashlib.md5(folder.encode()).hexdigest()[:6]}"

            else:
                 binner_name = base_name
                 # Check if this exact name already exists (from previous indexed fallback or manual input)
                 if any(name == binner_name for name in binner_mapping.values()):
                      binner_name = f"binner_{i+1:03d}" # Fallback to index

            binner_mapping[folder] = binner_name
            temp_name_counts[base_name] += 1 # Track base name usage


    # Save the mapping
    map_df = pd.DataFrame(list(binner_mapping.items()), columns=['FolderPath', 'BinnerName'])
    map_path = os.path.join(args.output_dir, "binner_mapping.tsv")
    map_df.to_csv(map_path, sep='\t', index=False)
    logger.info(f"Saved binner mapping to {map_path}")
    logger.info(f"Mapping used: {binner_mapping}")


    # Invert mapping for easy lookup: folder_path -> sanitized_name
    folder_to_sanitized_name = {folder: name for folder, name in binner_mapping.items()}


    # --- Generate Combinations (Size >= 2 only) ---
    logger.info(f"Generating combinations (size >= 2) from {len(args.bin_folders)} bin folders")
    combinations, _ = generate_bin_combinations(
        args.bin_folders,
        folder_to_sanitized_name, # Use the generated mapping
        args.max_combinations
    )
    if combinations:
        logger.info(f"Generated {len(combinations)} combinations (size >= 2)")
    else:
         logger.info("No combinations of size >= 2 generated.")


    # --- Run Refinement (binning_refiner.py) on combinations ---
    refined_bin_output_dirs = [] # Store paths to directories containing refined bins
    if not args.skip_refinement and combinations:
        logger.info("Starting refinement process on combinations")
        refined_bin_output_dirs = process_combinations_in_parallel(
            combinations,
            args.output_dir, # Pass base output dir
            args.binning_refiner_path,
            args.min_bin_size, # Pass min bin size for refiner (might differ from final QC)
            args.threads,
            args.batch_size,
            folder_to_sanitized_name, # Pass the mapping
            args.keep_temp # Pass keep_temp flag
        )
        logger.info(f"Refinement process completed. Found {len(refined_bin_output_dirs)} sets of refined bins from combinations.")
    elif args.skip_refinement:
        logger.info("Skipping refinement step (--skip-refinement).")
    elif not combinations:
        logger.info("Skipping refinement step as no combinations (size >= 2) were generated.")


    # --- Run CheckM2 (Unified Run on Originals + Refined) ---
    checkm_results = pd.DataFrame() # Initialize empty dataframe
    # Sources for CheckM2: original input folders + successful refinement output folders
    checkm_input_sources = args.bin_folders + refined_bin_output_dirs
    logger.info(f"Total sources for CheckM2: {len(checkm_input_sources)} ({len(args.bin_folders)} original + {len(refined_bin_output_dirs)} refined sets)")

    if not args.skip_checkm and checkm_input_sources:
        logger.info("Running unified CheckM2 on all relevant bin sets")
        checkm_results, bin_mapping_checkm = run_unified_checkm2(
             checkm_input_sources,
             output_dirs["checkm_results"], # Unified output dir
             args.checkm_db,
             args.threads,
             args.output_dir, # Pass base output dir for path logic
             folder_to_sanitized_name # Pass mapping for better source identification
        )
        if checkm_results.empty:
             logger.warning("CheckM2 run failed or produced no results. Bin selection will rely on size only.")
        else:
             logger.info("Unified CheckM2 run finished.")
    elif args.skip_checkm:
        logger.info("Skipping CheckM2 step (--skip-checkm).")
    else:
        logger.info("No bins found to run CheckM2 on (check input folders and refinement step).")


    # --- Consolidate Bins ---
    final_bin_count = 0 # Initialize count of final bins produced
    if not args.skip_consolidation:
        logger.info("Consolidating bins")

        # 1. Get contig mappings for ALL bins (originals + refined)
        logger.info("Collecting contig information from all sources...")
        all_bin_sources_for_consol = args.bin_folders + refined_bin_output_dirs
        contig_to_bins, contig_lengths, bin_info, all_bin_sets = get_contig_to_bin_mappings(
            all_bin_sources_for_consol,
            folder_to_sanitized_name # Pass the mapping
        )
        all_bin_ids = list(bin_info.keys())

        if not all_bin_ids:
            logger.error("No bins found for consolidation. Exiting.")
            # Clean up potentially empty dirs created if appropriate
            if not args.keep_temp:
                 temp_dir_to_remove = os.path.join(args.output_dir, "temp")
                 if os.path.isdir(temp_dir_to_remove): shutil.rmtree(temp_dir_to_remove)
            sys.exit(1)


        logger.info(f"Collected contig information for {len(all_bin_ids)} bins from {len(all_bin_sets)} sets.")

        # 2. Calculate Similarity
        logger.info("Calculating bin similarities")
        distance_matrix, similarity_matrix = calculate_bin_similarity(contig_to_bins, contig_lengths, all_bin_ids)

        # 3. Cluster Bins
        logger.info("Clustering bins")
        clusters, singletons = cluster_bins(distance_matrix, all_bin_ids, args.min_cluster_size, args.min_samples)

        if clusters is None:
            logger.error("Clustering failed. Cannot proceed with consolidation.")
            if not args.keep_temp:
                 temp_dir_to_remove = os.path.join(args.output_dir, "temp")
                 if os.path.isdir(temp_dir_to_remove): shutil.rmtree(temp_dir_to_remove)
            sys.exit(1)


        # Add singleton points back as singleton clusters
        singleton_clusters = {}
        singleton_label_start = 0
        numeric_labels = [k for k in clusters.keys() if isinstance(k, int)]
        if numeric_labels: singleton_label_start = max(numeric_labels) + 1
        for i, singleton_bin in enumerate(singletons):
             singleton_clusters[f"singleton_{singleton_label_start + i}"] = [singleton_bin] # Use string keys for singletons
        logger.info(f"Adding {len(singleton_clusters)} singleton points as singleton clusters.")
        clusters.update(singleton_clusters) # Add singleton clusters to main dict
        logger.info(f"Total potential clusters (including singletons): {len(clusters)}")


        # 4. Validate Clusters (Optional but recommended)
        logger.info("Validating clusters")
        valid_clusters, _ = validate_clusters(clusters, similarity_matrix, all_bin_ids, args.similarity_threshold)
        logger.info(f"Found {len(valid_clusters)} clusters after validation.")


        # 5. Visualize Clusters (Optional)
        if not args.skip_visualization and valid_clusters and len(all_bin_ids) > 1 : # UMAP needs > 1 sample
             logger.info("Creating bin clustering visualization")
             # Prepare data for visualization (e.g., using UMAP)
             try:
                  # Check for NaNs/Infs in distance matrix before UMAP
                  if np.isnan(distance_matrix).any() or np.isinf(distance_matrix).any():
                      logger.warning("Distance matrix contains NaN or Inf values. Skipping UMAP visualization.")
                  else:
                       # Reduce dimensionality using UMAP on the distance matrix
                       reducer = umap.UMAP(n_neighbors=max(2,min(15, len(all_bin_ids)-1)), # Adjust neighbors based on data size
                                           min_dist=0.1,
                                           metric='precomputed',
                                           random_state=42,
                                           n_components=2) # 2D for plotting
                       embedding = reducer.fit_transform(distance_matrix)

                       # Get cluster labels for all bins (map validated clusters back)
                       bin_labels_map = {}
                       for label, bins in valid_clusters.items():
                           for bin_id in bins:
                               bin_labels_map[bin_id] = str(label) # Ensure labels are strings for plotting

                       # Assign a label 'Split/Singleton' or similar if a bin ID isn't in any valid cluster
                       bin_labels = [bin_labels_map.get(b, 'Split/Singleton') for b in all_bin_ids]


                       plot_df = pd.DataFrame(embedding, columns=('UMAP_1', 'UMAP_2'))
                       plot_df['BinID'] = all_bin_ids
                       plot_df['ClusterLabel'] = bin_labels
                       plot_df['BinSet'] = [bin_info.get(b,{}).get('bin_set', 'Unknown') for b in all_bin_ids]
                       
                       # Static plot
                       try:
                           plt.figure(figsize=(12, 10)) # Keep figure size reasonable
                           unique_labels = sorted(list(set(bin_labels)))
                           n_labels = len(unique_labels)
                           logger.debug(f"Generating UMAP plot with {n_labels} unique cluster labels.")

                           # --- Improved Color Handling ---
                           # Use a qualitative colormap like tab20 or Paired, handle more colors if needed
                           if n_labels <= 10:
                               palette = sns.color_palette("tab10", n_colors=n_labels)
                           elif n_labels <= 20:
                               palette = sns.color_palette("tab20", n_colors=n_labels)
                           elif n_labels <= 24: # Use paired for up to 12 pairs
                               palette = sns.color_palette("Paired", n_colors=n_labels)
                           else: # Fallback for many labels - cycle 'tab20'
                               print(f"Warning: {n_labels} clusters exceed standard qualitative palettes. Colors will cycle.")
                               base_palette = sns.color_palette("tab20", n_colors=20)
                               palette = [base_palette[i % 20] for i in range(n_labels)]

                           # Ensure 'Split/Singleton' points (if present) are a consistent color (e.g., grey)
                           label_to_color = {label: palette[i] for i, label in enumerate(unique_labels) if label != 'Split/Singleton'}
                           if 'Split/Singleton' in unique_labels:
                               label_to_color['Split/Singleton'] = 'grey'
                               # Ensure grey is mapped correctly if 'Split/Singleton' wasn't the last label
                               ordered_unique_labels = [l for l in unique_labels if l != 'Split/Singleton'] + ['Split/Singleton']
                               final_palette = [label_to_color[l] for l in ordered_unique_labels]
                           else:
                               ordered_unique_labels = unique_labels
                               final_palette = palette

                           # Create scatter plot
                           scatter = sns.scatterplot(data=plot_df, x='UMAP_1', y='UMAP_2',
                                                     hue='ClusterLabel', hue_order=ordered_unique_labels,
                                                     style='BinSet', # Style by original binner/source
                                                     s=60,          # Point size
                                                     palette=final_palette, # Use the generated palette
                                                     legend='auto') # Let matplotlib handle legend initially

                           plt.title('Bin Clustering (HDBSCAN + UMAP)')
                           plt.grid(True, linestyle='--', alpha=0.6)

                           # --- Improved Legend Handling ---
                           current_legend = scatter.axes.get_legend()
                           if current_legend:
                                # Determine number of columns (aim for ~20-30 items per column max)
                                num_items = len(current_legend.get_texts())
                                ncol = max(1, int(np.ceil(num_items / 30)))

                                # Recreate legend outside plot area with smaller font and columns
                                scatter.legend(bbox_to_anchor=(1.02, 1), # Position right, top aligned
                                                loc='upper left',
                                                borderaxespad=0.,
                                                fontsize='xx-small', # Very small font
                                                ncol=ncol,          # Calculated columns
                                                title="Clusters & Sources", # Combined title
                                                title_fontsize='x-small')

                           cluster_plot_path = os.path.join(output_dirs["figures"], "bin_clustering.png")
                           plt.savefig(cluster_plot_path, dpi=150, bbox_inches='tight') # Use tight bbox
                           plt.close() # Close plot
                           logger.info(f"Saved bin clustering visualization to {cluster_plot_path}")

                       except Exception as e:
                           logger.warning(f"Could not create static UMAP plot: {e}")
                           # Attempt to close plot if it was opened
                           try: 
                               plt.close()
                           except: 
                               pass # Ignore errors during cleanup closing

                       # Interactive plot (optional, requires plotly)
                       try:
                           import plotly.express as px
                           fig = px.scatter(plot_df, x='UMAP_1', y='UMAP_2', color='ClusterLabel', symbol='BinSet', hover_data=['BinID'], title='Interactive Bin Clustering')
                           interactive_plot_path = os.path.join(output_dirs["figures"], "bin_clustering.html")
                           fig.write_html(interactive_plot_path)
                           logger.info(f"Saved interactive bin clustering visualization to {interactive_plot_path}")
                       except ImportError:
                           logger.info("Plotly not installed, skipping interactive visualization.")
                       except Exception as e:
                           logger.warning(f"Could not create interactive plot: {e}")

             except ImportError as e:
                  if 'umap' in str(e).lower():
                       logger.warning("UMAP-learn not installed. Skipping visualization. Install with: pip install umap-learn")
                  else:
                       logger.warning(f"Import error during visualization: {e}")
             except ValueError as e:
                  logger.warning(f"Value error during visualization (possibly UMAP issue): {e}")
             except Exception as e:
                  logger.warning(f"Could not create clustering visualization: {e}")


        # 6. Select Best Bin from each cluster and apply QC
        logger.info("Selecting best bins from each cluster and applying quality control...")
        best_bins_passing_qc = [] # List of best bin_ids that pass QC
        for label, cluster_members in valid_clusters.items():
             if not cluster_members: continue # Skip empty clusters (e.g., after validation)
             best_bin_id = select_best_bin(cluster_members, checkm_results, bin_info)
             if best_bin_id:
                  # Filter by quality thresholds AFTER selection
                  quality_ok = True
                  completeness = 'N/A'
                  contamination = 'N/A'
                  bin_length_kb = bin_info.get(best_bin_id,{}).get('length', 0) / 1000.0 # Length in kbp

                  # Check Length threshold first
                  if bin_length_kb < args.min_bin_size:
                       quality_ok = False
                       logger.debug(f"Bin {best_bin_id} (Cluster '{label}') failed length threshold ({bin_length_kb:.1f} kbp < {args.min_bin_size} kbp). Skipping.")
                  else:
                       # If length OK, check CheckM thresholds (if CheckM was run)
                       if not checkm_results.empty and 'original_path' in checkm_results.columns:
                             best_bin_path = bin_info.get(best_bin_id, {}).get('path')
                             if best_bin_path:
                                 checkm_row = checkm_results[checkm_results['original_path'] == best_bin_path]
                                 if not checkm_row.empty:
                                      comp = pd.to_numeric(checkm_row.iloc[0].get('Completeness'), errors='coerce')
                                      cont = pd.to_numeric(checkm_row.iloc[0].get('Contamination'), errors='coerce')
                                      completeness = f"{comp:.2f}" if pd.notna(comp) else 'N/A'
                                      contamination = f"{cont:.2f}" if pd.notna(cont) else 'N/A'
                                      # Check thresholds
                                      if pd.isna(comp) or pd.isna(cont) or comp < args.min_completion or cont > args.max_contamination:
                                          quality_ok = False
                                          logger.debug(f"Bin {best_bin_id} (Cluster '{label}') failed CheckM quality thresholds (C:{completeness}%, X:{contamination}%). Skipping.")
                                 else:
                                      logger.debug(f"No CheckM results for bin {best_bin_id} (Cluster '{label}'). Cannot apply CheckM QC.")
                             else:
                                  logger.warning(f"Cannot find path for best bin {best_bin_id} (Cluster '{label}') to check quality.")
                       elif not args.skip_checkm: # CheckM was supposed to run but results are empty
                            logger.warning(f"CheckM results are empty, cannot apply CheckM quality filter for bin {best_bin_id} (Cluster '{label}')")


                  if quality_ok:
                      logger.info(f"Bin {best_bin_id} (Cluster '{label}') passed all thresholds (Length: {bin_length_kb:.1f} kbp, C:{completeness}%, X:{contamination}%). Adding to final set.")
                      best_bins_passing_qc.append(best_bin_id)
                  # else: bin is skipped due to thresholds

             else:
                 logger.warning(f"No suitable bin found for cluster {label}.")

        logger.info(f"Selected {len(best_bins_passing_qc)} final bins meeting quality criteria.")

        # 7. Resolve Ambiguous Contigs among best bins passing QC
        logger.info("Resolving ambiguous contigs among selected bins")
        final_bins_contigs = resolve_ambiguous_contigs(
            best_bins_passing_qc, contig_to_bins, bin_info, args.keep_ambiguous, args.remove_ambiguous
        )
        logger.info(f"Number of bins remaining after ambiguity resolution: {len(final_bins_contigs)}")


        # 8. Write Final Bins
        logger.info("Writing final bins")
        final_bin_count = write_final_bins(
            final_bins_contigs,
            bin_info,
            output_dirs["final_bins"]
        )

        # Create final report using bins that passed QC and were potentially modified by ambiguity resolution
        logger.info("Creating final report")
        final_bin_ids_written = list(final_bins_contigs.keys()) # Use keys from the final dict
        create_final_report(
            final_bin_ids_written, # Report only on bins that were actually written
            checkm_results, # Pass the unified checkm results
            bin_info, # Pass original bin_info to lookup paths etc.
            os.path.join(args.output_dir, "final_bins_report.tsv")
        )

        # Plot bin quality using bins that passed QC and were written
        if not args.skip_visualization and final_bin_count > 0:
             logger.info("Creating final bin quality plot")
             plot_bin_quality(
                 final_bin_ids_written, # Plot only bins that were written
                 checkm_results, # Pass the unified checkm results
                 bin_info,
                 os.path.join(output_dirs["figures"], "final_bin_quality.png")
             )

    else: # Consolidation skipped
         logger.info("Skipping consolidation step (--skip-consolidation).")
         # If consolidation is skipped, 'final' bins are just the inputs + refined
         # No final report or specific QC plot is generated in this case.
         final_bin_count = 0


    # Run CheckM2 on final bins only if consolidation was performed and bins were produced
    if not args.skip_consolidation and not args.skip_checkm and final_bin_count > 0:
        logger.info("Running CheckM2 on final bins")
        run_checkm2_on_final_bins(
            output_dirs["final_bins"],
            os.path.join(args.output_dir, "final_bins_checkm2"), # Specific output dir for final checkm
            args.checkm_db,
            args.threads
        )
    elif args.skip_consolidation:
        logger.info("Consolidation skipped, not running CheckM2 on final bins.")
    elif final_bin_count == 0 and not args.skip_consolidation: # Consolidation ran but produced no bins
         logger.warning("No final bins were written, skipping final CheckM2 run.")


    # Clean up temporary files
    if not args.keep_temp:
        logger.info("Cleaning up temporary files")
        temp_dir_to_remove = os.path.join(args.output_dir, "temp")
        if os.path.isdir(temp_dir_to_remove):
             try:
                  shutil.rmtree(temp_dir_to_remove)
                  logger.info(f"Removed temporary directory: {temp_dir_to_remove}")
             except Exception as e:
                  logger.warning(f"Error cleaning up temporary directory {temp_dir_to_remove}: {str(e)}")
        else:
             logger.info("Temporary directory not found (or already removed), skipping cleanup.")

    logger.info("Enhanced bin refinement completed.")


if __name__ == "__main__":
    # Call main which handles arg parsing internally now
    main()