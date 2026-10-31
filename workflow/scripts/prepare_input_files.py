#!/usr/bin/env python

import os
import shutil
import subprocess
import sys

# --- Snakemake object stub for testing ---
# This section allows testing outside Snakemake
# Comment out or remove when running with Snakemake
# class SnakemakeStub:
#     def __init__(self):
#         # Create dummy input files for testing
#         self._create_dummy_files()
#         self.input = ["test_input.txt", "test_input.txt.gz", "test_input.txt.bz2"]
#         # Output names MUST ideally end in .gz as defined in the Snakefile rule
#         self.output = ["output.txt.gz", "output.gz.gz", "output.bz2.gz"]
#         self.threads = 4
#
#     def _create_dummy_files(self):
#         print("Creating dummy files for testing...")
#         with open("test_input.txt", "w") as f:
#             f.write("This is a test file.\n" * 10)
#         subprocess.run("gzip -k test_input.txt", shell=True, check=True)
#         subprocess.run("bzip2 -k test_input.txt", shell=True, check=True)
#         print("Dummy files created: test_input.txt, test_input.txt.gz, test_input.txt.bz2")
#
# # Comment out the next line when running with Snakemake
# snakemake = SnakemakeStub()
# --- End Snakemake stub ---


# Validate inputs and outputs
if not hasattr(snakemake, 'input') or not hasattr(snakemake, 'output') or not hasattr(snakemake, 'threads'):
     print("Error: Script expects snakemake object with 'input', 'output', and 'threads' attributes.", file=sys.stderr)
     sys.exit(1)

if len(snakemake.input) == 0:
    print("Warning: No input files provided.", file=sys.stderr)
    # Allow script to finish successfully if no inputs, consistent with Snakemake's handling of empty inputs.
    sys.exit(0)

if len(snakemake.input) != len(snakemake.output):
    raise ValueError(f"Input and output file lists have different lengths: "
                     f"{len(snakemake.input)} inputs ({', '.join(snakemake.input)}) vs "
                     f"{len(snakemake.output)} outputs ({', '.join(snakemake.output)})")

# Ensure all designated output files end with .gz or .gzip
# This is mainly a sanity check; the Snakemake rule should define this.
for ofile in snakemake.output:
    if not ofile.lower().endswith((".gz", ".gzip")):
        print(f"Warning: Output file '{ofile}' does not end with .gz or .gzip. "
              f"The script will still create it as a gzipped file.", file=sys.stderr)

# Process each input/output pair
for ifile, ofile in zip(snakemake.input, snakemake.output):
    print(f"Processing: {ifile} => {ofile} (Ensuring output is gzipped)")

    # Get the lowercased extension of the input file
    try:
        ext = os.path.splitext(ifile)[-1].lower()
    except IndexError:
        ext = "" # Handle files with no extension

    try:
        # Case 1: Input is already gzipped (.gz, .gzip)
        if ext in ['.gz', '.gzip']:
            print(f"  Input '{ifile}' is gzipped. Copying directly to '{ofile}'.")
            # If input is already correctly gzipped, just copy it. Most efficient.
            shutil.copy(ifile, ofile)

        # Case 2: Input is bzip2 compressed (.bz2, .bzip2)
        elif ext in ['.bz2', '.bzip2']:
            print(f"  Input '{ifile}' is bzip2. Decompressing and re-compressing to gzip '{ofile}'.")
            # Decompress bzip2 using pbzip2 and pipe to pigz for gzip compression
            # Use bash -c '...' to handle the pipe correctly within shell=True
            # Using f-string for clarity, ensure proper quoting if filenames have spaces (though Snakemake usually handles this)
            cmd = f"pbzip2 --decompress --stdout -p{snakemake.threads} '{ifile}' | pigz -p{snakemake.threads} > '{ofile}'"
            print(f"  Running command: {cmd}")
            # Using subprocess.run for better control and error handling
            # check=True raises CalledProcessError on non-zero exit status
            # shell=True is needed for the pipe '|'
            # executable='/bin/bash' is good practice when using shell=True with complex commands/pipes
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')

        # Case 3: Input is uncompressed (or unknown compression type)
        else:
            print(f"  Input '{ifile}' appears uncompressed. Compressing with pigz to '{ofile}'.")
            # Compress the input file using pigz, writing to stdout (-c) and redirecting
            cmd = f"pigz -p{snakemake.threads} -c '{ifile}' > '{ofile}'"
            print(f"  Running command: {cmd}")
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')

    except subprocess.CalledProcessError as e:
        print(f"\nError processing {ifile} to {ofile}.", file=sys.stderr)
        print(f"Command failed with exit code {e.returncode}:", file=sys.stderr)
        print(f"Command: {e.cmd}", file=sys.stderr)
        # Attempt to remove potentially incomplete/corrupt output file
        if os.path.exists(ofile):
            print(f"Cleaning up incomplete output file: {ofile}", file=sys.stderr)
            try:
                os.remove(ofile)
            except OSError as rm_err:
                print(f"Warning: Failed to remove incomplete output file '{ofile}': {rm_err}", file=sys.stderr)
        sys.exit(1) # Exit script with error status

    except Exception as e:
        print(f"\nAn unexpected error occurred processing {ifile} to {ofile}: {e}", file=sys.stderr)
        # Attempt cleanup here too
        if os.path.exists(ofile):
            print(f"Cleaning up potentially incomplete output file: {ofile}", file=sys.stderr)
            try:
                os.remove(ofile)
            except OSError as rm_err:
                 print(f"Warning: Failed to remove incomplete output file '{ofile}': {rm_err}", file=sys.stderr)
        sys.exit(1) # Exit script with error status


print("\nScript finished successfully.")

# --- Cleanup for testing stub ---
# Comment out or remove when running with Snakemake
# if isinstance(snakemake, SnakemakeStub):
#    print("\nCleaning up dummy and output files...")
#    files_to_remove = ["test_input.txt", "test_input.txt.gz", "test_input.txt.bz2"] + snakemake.output
#    for f in files_to_remove:
#        if os.path.exists(f):
#            os.remove(f)
#            print(f"Removed: {f}")
# --- End cleanup ---